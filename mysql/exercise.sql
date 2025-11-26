/*
 * SQL Practice Exercises
 * Comprehensive set of SQL queries covering various operations
 * Database: sakila (MySQL sample database), employees (sample database)
 */

-- ============================================
-- SET 1: BASIC SELECT OPERATIONS
-- ============================================

-- 1a. Select all columns from the actor table
SELECT * FROM actor;

-- 1b. Select only the last_name column from the actor table
SELECT last_name FROM actor;

-- ============================================
-- SET 2: DISTINCT VALUES
-- ============================================

-- 2a. Select all distinct last names from the actor table
SELECT DISTINCT last_name FROM actor;

-- 2b. Select all distinct postal codes from the address table
SELECT DISTINCT postal_code FROM address;

-- 2c. Select all distinct ratings from the film table
SELECT DISTINCT rating FROM film;

-- ============================================
-- SET 3: WHERE CLAUSE FILTERING
-- ============================================

-- 3a. Select films with length >= 180 minutes (3 hours or longer)
SELECT title, description, rating, length
FROM film
WHERE length >= 180;

-- 3b. Select payments made on or after 05/27/2005
SELECT payment_id, amount, payment_date
FROM payment
WHERE payment_date >= '2005-05-27 00:00:00';

-- 3c. Select payments made specifically on 05/27/2005
SELECT payment_id, amount, payment_date
FROM payment
WHERE payment_date >= '2005-05-27 00:00:00' AND payment_date < '2005-05-28 00:00:00';

-- 3d. Select customers where last name begins with "S" and first name ends with "N"
SELECT *
FROM customer
WHERE last_name LIKE 'S%' AND first_name LIKE '%N';

-- 3e. Select categories where ID > 4 and name begins with C, S, or T
SELECT *
FROM category
WHERE category_id > 4
  AND (name LIKE 'C%' OR name LIKE 'S%' OR name LIKE 'T%');

-- 3f. Select all staff columns except password where password exists
SELECT staff_id, first_name, last_name, address_id, picture, email, store_id, active, username, last_update
FROM staff
WHERE password IS NOT NULL;

-- ============================================
-- SET 4: SORTING AND ORDERING
-- ============================================

-- 4a. Select all films ordered by length in ascending order
SELECT *
FROM film
ORDER BY length ASC;

-- 4b. Select distinct ratings ordered in descending order
SELECT DISTINCT rating
FROM film
ORDER BY rating DESC;

-- 4c. Select top 20 payments ordered by amount (descending)
SELECT payment_date, amount
FROM payment
ORDER BY amount DESC
LIMIT 20;

-- ============================================
-- SET 5: PATTERN MATCHING WITH LIKE
-- ============================================

-- 5a. Select films where description begins with "A Thoughtful"
SELECT title, description, release_year
FROM film
WHERE description LIKE 'A Thoughtful%';

-- 5b. Select films where description ends with "Boat"
SELECT title, description, rental_duration
FROM film
WHERE description LIKE '%Boat';

-- 5c. Select films where description contains "Database" and length > 180 minutes
SELECT title, length, description, rental_rate
FROM film
WHERE description LIKE '%Database%' AND length > 180;

-- ============================================
-- SET 6: IN OPERATOR
-- ============================================

-- 6a. Select addresses in specific districts
SELECT phone, district
FROM address
WHERE district IN ('California', 'England', 'Taipei', 'West Java');

-- 6b. Select payments made on specific dates
SELECT payment_id, amount, payment_date
FROM payment
WHERE (payment_date >= '2005-05-25 00:00:00' AND payment_date < '2005-05-26 00:00:00')
   OR (payment_date >= '2005-05-27 00:00:00' AND payment_date < '2005-05-28 00:00:00')
   OR (payment_date >= '2005-05-29 00:00:00' AND payment_date < '2005-05-30 00:00:00');

-- 6c. Select films with specific ratings
SELECT *
FROM film
WHERE rating IN ('G', 'PG-13', 'NC-17');

-- ============================================
-- SET 7: BETWEEN OPERATOR
-- ============================================

-- 7a. Retrieve customers with customer_id between 10 and 50
SELECT *
FROM customer
WHERE customer_id BETWEEN 10 AND 50;

-- 7b. Find films with length between 90 and 120 minutes
SELECT *
FROM film
WHERE length BETWEEN 90 AND 120;

-- 7c. Get payments with amounts between 5.00 and 7.00
SELECT *
FROM payment
WHERE amount BETWEEN 5.00 AND 7.00;

-- 7d. List films with release_year between 2005 and 2010
SELECT *
FROM film
WHERE release_year BETWEEN 2005 AND 2010;

-- ============================================
-- SET 8: ADVANCED FILTERING
-- ============================================

-- 2.1. Find rentals from 2005 excluding May
SELECT *
FROM rental
WHERE rental_date >= '2005-01-01 00:00:00' AND rental_date < '2006-01-01 00:00:00'
  AND (rental_date < '2005-05-01 00:00:00' OR rental_date >= '2005-06-01 00:00:00');

-- 2.2. Find movies not rated G or PG
SELECT *
FROM film
WHERE rating NOT IN ('G', 'PG');

-- ============================================
-- SET 9: GROUP BY AND AGGREGATION
-- ============================================

-- 2.3. Count total rentals per staff member
SELECT staff_id, COUNT(*) AS total_rentals
FROM rental
GROUP BY staff_id;

-- 2.4. Count distinct films per store
SELECT store_id, COUNT(DISTINCT film_id) AS total_films
FROM inventory
GROUP BY store_id;

-- 2.5. Find rental dates with more than 100 rentals
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
HAVING COUNT(*) > 100;

-- 2.6. Find actors who appeared in more than 20 films
SELECT actor_id, COUNT(*) AS film_count
FROM film_actor
GROUP BY actor_id
HAVING COUNT(*) > 20;

-- 2.7. Find customers who spent more than $200
SELECT customer_id, SUM(amount) AS total_spent
FROM payment
GROUP BY customer_id
HAVING SUM(amount) > 200;

-- 2.8. Convert payment amounts from USD to BDT (1 USD = 120 BDT)
SELECT customer_id,
       SUM(amount) AS total_usd,
       SUM(amount) * 120 AS total_bdt
FROM payment
GROUP BY customer_id;

-- ============================================
-- SET 10: MULTI-COLUMN SORTING
-- ============================================

-- 2.9. List R and NC-17 rated films sorted alphabetically
SELECT title, rating
FROM film
WHERE rating IN ('R', 'NC-17')
ORDER BY title ASC;

-- 2.10. Sort films by rental_rate (DESC) then length (ASC)
SELECT title, rental_rate, length
FROM film
ORDER BY rental_rate DESC, length ASC;

-- 2.11. Get first 10 films sorted alphabetically by title
SELECT title, description, release_year
FROM film
ORDER BY title ASC
LIMIT 10;

-- 2.12. Find inventory items rented less than 2 times
SELECT inventory_id
FROM rental
GROUP BY inventory_id
HAVING COUNT(*) < 2;

-- ============================================
-- SET 11: JOINS
-- ============================================

-- 4.1. Find customers who made more than 5 payments
SELECT customer_id
FROM payment
GROUP BY customer_id
HAVING COUNT(*) > 5;

-- 4.2. Find rental dates with more than 100 rentals
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
HAVING COUNT(*) > 100;

-- 4.3. Count total payments per staff member
SELECT staff_id, COUNT(*) AS total_payments
FROM payment
GROUP BY staff_id;

-- 4.4. Find top 5 busiest rental days
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
ORDER BY COUNT(*) DESC
LIMIT 5;

-- 4.5. Find maximum replacement cost per rating
SELECT rating, MAX(replacement_cost) AS max_replacement_cost
FROM film
GROUP BY rating;

-- 4.6. Find films with rental_rate greater than 3
SELECT title, rental_rate
FROM film
WHERE rental_rate > 3;

-- 4.7. Count rentals by year and month
SELECT YEAR(rental_date) AS year, MONTH(rental_date) AS month, COUNT(*) AS total_rentals
FROM rental
GROUP BY YEAR(rental_date), MONTH(rental_date)
ORDER BY year, month;

-- 4.8. Find customers who made payments on 2005-07-30
SELECT DISTINCT c.first_name, c.last_name
FROM customer c
JOIN payment p ON p.customer_id = c.customer_id
WHERE DATE(p.payment_date) = '2005-07-30';

-- 4.9. Find postal code for customer ID 2
SELECT a.postal_code
FROM customer c
JOIN address a ON a.address_id = c.address_id
WHERE c.customer_id = 2;

-- 4.10. Find customers who paid more than $50 in total
SELECT c.first_name, c.last_name, SUM(p.amount) AS total_paid
FROM customer c
JOIN payment p ON p.customer_id = c.customer_id
GROUP BY c.customer_id, c.first_name, c.last_name
HAVING SUM(p.amount) > 50;

-- ============================================
-- SET 12: EMPLOYEES SCHEMA QUERIES
-- Note: These queries require the 'employees' sample database
-- ============================================

-- 5.1. Find employees hired after January 1, 2000
SELECT first_name, last_name
FROM employees
WHERE hire_date > '2000-01-01';

-- 5.2. Find employees whose last name starts with 'M'
SELECT first_name, last_name
FROM employees
WHERE last_name LIKE 'M%';

-- 5.3. Get top 10 highest-paid employees
SELECT employee_id
FROM employees
ORDER BY salary DESC
LIMIT 10;

-- 5.4. Count employees by gender
SELECT gender, COUNT(*) AS total_employees
FROM employees
GROUP BY gender;

-- 5.5. Count employees hired per year
SELECT YEAR(hire_date) AS hire_year, COUNT(*) AS hires
FROM employees
GROUP BY YEAR(hire_date);

-- 5.6. Display employee ID, first name, and last name
SELECT employee_id, first_name, last_name
FROM employees;

-- 5.7. Display employee ID, job title, and salary
SELECT employee_id, job_title, salary
FROM employees;

-- 5.8. Find average salary per job title
SELECT job_title, AVG(salary) AS avg_salary
FROM employees
GROUP BY job_title;
