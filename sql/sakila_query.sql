-- --- Exercise Set 1: Fundamental Queries (Assuming Sakila-like schema) ---

-- 1a. Select all columns from the actor table.
SELECT * FROM actor;

-- 1b. Select only the last_name column from the actor table.
SELECT last_name FROM actor;

-- 2a. Select all distinct last names from the actor table.
SELECT DISTINCT last_name FROM actor;

-- 2b. Select all distinct postal codes from the address table.
SELECT DISTINCT postal_code FROM address;

-- 2c. Select all distinct ratings from the film table.
SELECT DISTINCT rating FROM film;

-- 3a. Select the title, description, rating, movie length columns from the films table that last 3 hours or longer (180 minutes).
SELECT title, description, rating, length FROM film WHERE length >= 180;

-- 3b. Select the payment id, amount, and payment date columns from the payments table for payments made on or after 05/27/2005.
SELECT payment_id, amount, payment_date FROM payment WHERE payment_date >= '2005-05-27';

-- 3c. Select the primary key, amount, and payment date columns from the payment table for payments made on 05/27/2005.
SELECT payment_id, amount, payment_date FROM payment WHERE DATE(payment_date) = '2005-05-27';

-- 3d. Select all columns from the customer table for rows where the last name begins with "S" and the first name ends with "N".
SELECT * FROM customer WHERE last_name LIKE 'S%' AND first_name LIKE '%N';

-- 3e. Select all columns from the category table for rows where the primary key is greater than 4 and the name field begins with "C", "S", or "T".
SELECT * FROM category WHERE category_id > 4 AND name REGEXP '^[CST]'; -- Using REGEXP for common SQL dialects
-- Alternative using standard SQL:
-- SELECT * FROM category WHERE category_id > 4 AND (name LIKE 'C%' OR name LIKE 'S%' OR name LIKE 'T%');

-- 3f. Select all columns minus the password column from the staff table for rows that contain a password (assuming the password column is named 'password').
SELECT staff_id, first_name, last_name, address_id, picture, email, store_id, active, username, last_update FROM staff WHERE password IS NOT NULL AND password <> '';

-- 4a. Select all columns from the film table and order rows by the length field in ascending order.
SELECT * FROM film ORDER BY length ASC;

-- 4b. Select all distinct ratings from the film table ordered by rating in descending order.
SELECT DISTINCT rating FROM film ORDER BY rating DESC;

-- 4c. Select the payment date and amount columns from the payment table for the first 20 payments ordered by payment amount in descending order.
SELECT payment_date, amount FROM payment ORDER BY amount DESC LIMIT 20;

-- 5a. Select the following columns from the film table for rows where the description begins with "A Thoughtful".
SELECT title, description, release_year FROM film WHERE description LIKE 'A Thoughtful%';

-- 5b. Select the following columns from the film table for rows where the description ends with the word "Boat".
SELECT title, description, rental_duration FROM film WHERE description LIKE '%Boat';

-- 5c. Select the following columns from the film table where the description contains the word "Database" and the length of the film is greater than 3 hours (180 minutes).
SELECT title, length, description, rental_rate FROM film WHERE description LIKE '%Database%' AND length > 180;

-- 6a. Select the phone and district columns from the address table for addresses in California, England, Taipei, or West Java.
SELECT phone, district FROM address WHERE district IN ('California', 'England', 'Taipei', 'West Java');

-- 6b. Select the payment id, amount, and payment date columns from the payment table for payments made on 05/25/2005, 05/27/2005, and 05/29/2005.
SELECT payment_id, amount, payment_date FROM payment WHERE DATE(payment_date) IN ('2005-05-25', '2005-05-27', '2005-05-29');

-- 6c. Select all columns from the film table for films rated G, PG-13, or NC-17.
SELECT * FROM film WHERE rating IN ('G', 'PG-13', 'NC-17');

-- 7a. Retrieve all customers with customer_id between 10 and 50.
SELECT * FROM customer WHERE customer_id BETWEEN 10 AND 50;

-- 7b. Find all films with length between 90 and 120 minutes.
SELECT * FROM film WHERE length BETWEEN 90 AND 120;

-- 7c. Get payments with amounts between 5.00 and 7.00.
SELECT * FROM payment WHERE amount BETWEEN 5.00 AND 7.00;

-- 7d. List all films with release_year between 2005 and 2010.
SELECT * FROM film WHERE release_year BETWEEN 2005 AND 2010;

-- --- Exercise Set 2: Joins and Aggregation (Assuming Sakila-like schema) ---

-- 1. Find all rentals made in 2005 but exclude those made in May.
SELECT * FROM rental
WHERE YEAR(rental_date) = 2005 AND MONTH(rental_date) <> 5;

-- 2. Find all movies that are neither rated G nor PG.
SELECT * FROM film
WHERE rating NOT IN ('G', 'PG');

-- 3. Find the staff ID and the total number of rentals processed by each staff member.
SELECT staff_id, COUNT(rental_id) AS total_rentals
FROM rental
GROUP BY staff_id;

-- 4. Find the store ID and the total number of films available in each store.
SELECT inventory.store_id, COUNT(inventory.inventory_id) AS total_films_available
FROM inventory
GROUP BY inventory.store_id;

-- 5. Find rental dates with more than 100 rentals.
SELECT DATE(rental_date) AS rental_day, COUNT(rental_id) AS total_rentals
FROM rental
GROUP BY rental_day
HAVING total_rentals > 100;

-- 6. Find actors who have acted in more than 20 films.
SELECT fa.actor_id, a.first_name, a.last_name, COUNT(fa.film_id) AS film_count
FROM film_actor fa
JOIN actor a ON fa.actor_id = a.actor_id
GROUP BY fa.actor_id, a.first_name, a.last_name
HAVING film_count > 20;

-- 7. Find customers who spent more than $200 on rentals.
SELECT p.customer_id, c.first_name, c.last_name, SUM(p.amount) AS total_spent
FROM payment p
JOIN customer c ON p.customer_id = c.customer_id
GROUP BY p.customer_id, c.first_name, c.last_name
HAVING total_spent > 200;

-- 8. Convert the total payment amount for each customer from USD to BDT using a conversion rate of 120 BDT per USD.
SELECT customer_id,
       SUM(amount) AS total_payment_usd,
       SUM(amount) * 120 AS total_payment_bdt
FROM payment
GROUP BY customer_id;

-- 9. Retrieve a list of films that are rated "R" or "NC-17". Display the results in two columns: title and rating, sorted alphabetically by title.
SELECT title, rating
FROM film
WHERE rating IN ('R', 'NC-17')
ORDER BY title ASC;

-- 10. Retrieve the title, rental_rate, and length of films. Sort the results in descending order of rental_rate, and then by length in ascending order.
SELECT title, rental_rate, length
FROM film
ORDER BY rental_rate DESC, length ASC;

-- 11. Retrieve the title, description, and release_year of films. Sort the results alphabetically by title, and display only the first 10 films.
SELECT title, description, release_year
FROM film
ORDER BY title ASC
LIMIT 10;

-- 12. Find inventory ids rented less than 2 times.
SELECT inventory_id
FROM rental
GROUP BY inventory_id
HAVING COUNT(rental_id) < 2;

-- --- Exercise Set 3: Employees Database Queries (Using `employees` schema tables like employees, titles, salaries) ---

-- 1. Display the first name, last name, and hire date of all employees.
SELECT first_name, last_name, hire_date FROM employees;

-- 2. List all distinct job titles available.
SELECT DISTINCT title FROM titles;

-- 3. Find the first and last names of employees who were hired after January 1, 2000.
SELECT first_name, last_name FROM employees WHERE hire_date > '2000-01-01';

-- 4. Display the first and last names of employees whose last name starts with 'M'.
SELECT first_name, last_name FROM employees WHERE last_name LIKE 'M%';

-- 5. Retrieve the employee IDs of the top 10 highest-paid employees. Sort the results in descending order of salary.
SELECT emp_no FROM salaries ORDER BY salary DESC LIMIT 10;

-- 6. Find the total number of employees grouped by gender.
SELECT gender, COUNT(emp_no) AS total_employees FROM employees GROUP BY gender;

-- 7. Count how many employees were hired in each year.
SELECT YEAR(hire_date) AS hire_year, COUNT(emp_no) AS total_hired
FROM employees
GROUP BY hire_year
ORDER BY hire_year;

-- 8. Write a query to display the employee ID, first name, and last name for all employees.
SELECT emp_no, first_name, last_name FROM employees;

-- 9. Write a query to display the employee ID, job title, and salary of all employees (assuming current titles and salaries).
SELECT
    e.emp_no,
    t.title AS job_title,
    s.salary
FROM employees e
JOIN titles t ON e.emp_no = t.emp_no
JOIN salaries s ON e.emp_no = s.emp_no
WHERE t.to_date = '9999-01-01' -- Filter for current title
AND s.to_date = '9999-01-01'; -- Filter for current salary

-- 10. Find the average salary of employees grouped by their current job title.
SELECT
    t.title AS job_title,
    AVG(s.salary) AS average_salary
FROM salaries s
JOIN titles t ON s.emp_no = t.emp_no
WHERE s.to_date = '9999-01-01' -- Filter for current salary
AND t.to_date = '9999-01-01' -- Filter for current title
GROUP BY t.title;