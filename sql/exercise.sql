-- 1a
SELECT * FROM actor;

-- 1b
SELECT last_name FROM actor;

-- 2a
SELECT DISTINCT last_name FROM actor;

-- 2b
SELECT DISTINCT postal_code FROM address;

-- 2c
SELECT DISTINCT rating FROM film;

-- 3a
SELECT title, description, rating, length
FROM film
WHERE length >= 180;

-- 3b
SELECT payment_id, amount, payment_date
FROM payment
WHERE payment_date >= '2005-05-27';

-- 3c
SELECT payment_id, amount, payment_date
FROM payment
WHERE DATE(payment_date) = '2005-05-27';

-- 3d
SELECT *
FROM customer
WHERE last_name LIKE 'S%' AND first_name LIKE '%N';

-- 3e
SELECT *
FROM category
WHERE category_id > 4
  AND (name LIKE 'C%' OR name LIKE 'S%' OR name LIKE 'T%');

-- 3f
SELECT staff_id, first_name, last_name, address_id, picture, email, store_id, active, username, last_update
FROM staff
WHERE password IS NOT NULL;

-- 4a
SELECT *
FROM film
ORDER BY length ASC;

-- 4b
SELECT DISTINCT rating
FROM film
ORDER BY rating DESC;

-- 4c
SELECT payment_date, amount
FROM payment
ORDER BY amount DESC
LIMIT 20;

-- 5a
SELECT title, description, release_year
FROM film
WHERE description LIKE 'A Thoughtful%';

-- 5b
SELECT title, description, rental_duration
FROM film
WHERE description LIKE '%Boat';

-- 5c
SELECT title, length, description, rental_rate
FROM film
WHERE description LIKE '%Database%' AND length > 180;

-- 6a
SELECT phone, district
FROM address
WHERE district IN ('California', 'England', 'Taipei', 'West Java');

-- 6b
SELECT payment_id, amount, payment_date
FROM payment
WHERE DATE(payment_date) IN ('2005-05-25', '2005-05-27', '2005-05-29');

-- 6c
SELECT *
FROM film
WHERE rating IN ('G', 'PG-13', 'NC-17');

-- 7a
SELECT *
FROM customer
WHERE customer_id BETWEEN 10 AND 50;

-- 7b
SELECT *
FROM film
WHERE length BETWEEN 90 AND 120;

-- 7c
SELECT *
FROM payment
WHERE amount BETWEEN 5.00 AND 7.00;

-- 7d
SELECT *
FROM film
WHERE release_year BETWEEN 2005 AND 2010;

-- 2.1
SELECT *
FROM rental
WHERE YEAR(rental_date) = 2005
  AND MONTH(rental_date) <> 5;

-- 2.2
SELECT *
FROM film
WHERE rating NOT IN ('G', 'PG');

-- 2.3
SELECT staff_id, COUNT(*) AS total_rentals
FROM rental
GROUP BY staff_id;

-- 2.4
SELECT store_id, COUNT(DISTINCT film_id) AS total_films
FROM inventory
GROUP BY store_id;

-- 2.5
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
HAVING COUNT(*) > 100;

-- 2.6
SELECT actor_id, COUNT(*) AS film_count
FROM film_actor
GROUP BY actor_id
HAVING COUNT(*) > 20;

-- 2.7
SELECT customer_id, SUM(amount) AS total_spent
FROM payment
GROUP BY customer_id
HAVING SUM(amount) > 200;

-- 2.8
SELECT customer_id,
       SUM(amount) AS total_usd,
       SUM(amount) * 120 AS total_bdt
FROM payment
GROUP BY customer_id;

-- 2.9
SELECT title, rating
FROM film
WHERE rating IN ('R', 'NC-17')
ORDER BY title ASC;

-- 2.10
SELECT title, rental_rate, length
FROM film
ORDER BY rental_rate DESC, length ASC;

-- 2.11
SELECT title, description, release_year
FROM film
ORDER BY title ASC
LIMIT 10;

-- 2.12
SELECT inventory_id
FROM rental
GROUP BY inventory_id
HAVING COUNT(*) < 2;

-- 3.1
SELECT first_name, last_name, hire_date
FROM employees;

-- 3.2
SELECT DISTINCT job_title
FROM employees;

-- 3.3
SELECT first_name, last_name
FROM employees
WHERE hire_date > '2000-01-01';

-- 3.4
SELECT first_name, last_name
FROM employees
WHERE last_name LIKE 'M%';

-- 3.5
SELECT employee_id
FROM employees
ORDER BY salary DESC
LIMIT 10;

-- 3.6
SELECT gender, COUNT(*) AS total_employees
FROM employees
GROUP BY gender;

-- 3.7
SELECT YEAR(hire_date) AS hire_year, COUNT(*) AS hires
FROM employees
GROUP BY YEAR(hire_date);

-- 3.8
SELECT employee_id, first_name, last_name
FROM employees;

-- 3.9
SELECT employee_id, job_title, salary
FROM employees;

-- 3.10
SELECT job_title, AVG(salary) AS avg_salary
FROM employees
GROUP BY job_title;

-- 4.1
SELECT customer_id
FROM payment
GROUP BY customer_id
HAVING COUNT(*) > 5;

-- 4.2
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
HAVING COUNT(*) > 100;

-- 4.3
SELECT staff_id, COUNT(*) AS total_payments
FROM payment
GROUP BY staff_id;

-- 4.4
SELECT DATE(rental_date) AS rental_day, COUNT(*) AS total_rentals
FROM rental
GROUP BY DATE(rental_date)
ORDER BY COUNT(*) DESC
LIMIT 5;

-- 4.5
SELECT rating, MAX(replacement_cost) AS max_replacement_cost
FROM film
GROUP BY rating;

-- 4.6
SELECT title, rental_rate
FROM film
WHERE rental_rate > 3;

-- 4.7
SELECT YEAR(rental_date) AS year, MONTH(rental_date) AS month, COUNT(*) AS total_rentals
FROM rental
GROUP BY YEAR(rental_date), MONTH(rental_date)
ORDER BY year, month;

-- 4.8
SELECT DISTINCT c.first_name, c.last_name
FROM customer c
JOIN payment p ON p.customer_id = c.customer_id
WHERE DATE(p.payment_date) = '2005-07-30';

-- 4.9
SELECT a.postal_code
FROM customer c
JOIN address a ON a.address_id = c.address_id
WHERE c.customer_id = 2;

-- 4.10
SELECT c.first_name, c.last_name, SUM(p.amount) AS total_paid
FROM customer c
JOIN  payment p ON p.customer_id = c.customer_id
GROUP BY c.customer_id, c.first_name, c.last_name
HAVING SUM(p.amount) > 50;


-- employees SCHEMA

-- 5.1
SELECT first_name, last_name
FROM employees
WHERE hire_date > '2000-01-01';

-- 5.2
SELECT first_name, last_name
FROM employees
WHERE last_name LIKE 'M%';

-- 5.3
SELECT employee_id
FROM employees
ORDER BY salary DESC
LIMIT 10;

-- 5.4
SELECT gender, COUNT(*) AS total_employees
FROM employees
GROUP BY gender;

-- 5.5
SELECT YEAR(hire_date) AS hire_year, COUNT(*) AS hires
FROM employees
GROUP BY YEAR(hire_date);

-- 5.6
SELECT employee_id, first_name, last_name
FROM employees;

-- 5.7
SELECT employee_id, job_title, salary
FROM employees;

-- 5.8
SELECT job_title, AVG(salary) AS avg_salary
FROM employees
GROUP BY job_title;
