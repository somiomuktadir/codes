-- SQL Table

CREATE SCHEMA employee;
USE employee;

-- Creating a table
CREATE TABLE person(
	person_id SMALLINT UNSIGNED PRIMARY KEY AUTO_INCREMENT,
    fname VARCHAR(20),
    lname VARCHAR(20),
    eye_color ENUM("BR","BL", "GR"),
    birth_date DATE,
    address VARCHAR(20),
    email_info VARCHAR(20)
);

-- Inserting information
INSERT INTO person
(person_id, fname, lname, eye_color, birth_date, address, email_info) VALUES

("1", "William", "Turner", "BR", "2005-09-14", "MeowMeow Street", "wturner@gmail.com"),
(NULL, "Max", "Verstappen", "GR", "2005-09-15", "Red Bull Street", "maxmaxmax@gmail.com");


-- Finding information
SELECT * FROM person;

-- Finding information with a condition
SELECT * FROM person
WHERE person_id = 1;

-- Updating Data
UPDATE person
SET address="Meow Meow Street"
WHERE person_id=1;




-- SQL Query

-- SELECT 
SELECT name AS language_name, last_update AS updated_at
FROM sakila.language;

-- WHERE
SELECT film_id, title FROM film WHERE
(rating ="G" AND rental_duration >=7) OR (rating="NC-17" AND rental_duration >=4);

-- ORDER
SELECT customer_id, rental_date FROM rental
WHERE DATE(rental_date) = "2005-05-24"
ORDER BY rental_date DESC;

-- Range
SELECT customer_id, rental_date FROM rental
WHERE DATE(rental_date) <= "2005-06-15"
	AND DATE(rental_date) <= "2005-06-14";
    
SELECT customer_id, payment_date, amount FROM payment
WHERE amount BETWEEN 10.0 AND 11.99;

SELECT title, description FROM film
WHERE rating NOT IN ("PG-13", "G");


-- LIMIT & LIKE
SELECT address FROM address
WHERE address LIKE "%Street%"
LIMIT 10;

-- Basic Functions
-- SUM, MAX, COUNT, AVG

SELECT DATE(payment_date) AS paydate, SUM(amount) AS total_payment
FROM payment;

-- Grouping and Filtering
SELECT customer_id, COUNT(*) AS rentals
FROM rental
GROUP BY customer_id
HAVING rentals >= 40;





-- Subqueries


-- Find customer's name who made payments greater than 11
SELECT first_name, last_name FROM customer
WHERE customer_id IN (
	SELECT customer_id
    FROM sakila.payment
	WHERE amount > 11
);

-- Find all films that were never rented
SELECT title FROM film
WHERE film_id IN (
		SELECT film_id FROM inventory
        WHERE inventory_id NOT IN (
			SELECT inventory_id FROM rental
	)
);


-- Inner Join

-- Find address if customer id 1
SELECT address FROM address
INNER JOIN customer ON address.address_id=customer.address_id
WHERE customer.customer_id = 1;

-- LEFT JOIN and RIGHT JOIN

-- LEFT JOIN : table on the left side determines number of rows in the result set
-- RIGHT JOIN : table on the right side determines number of rows in the result set

SELECT film.film_id, film.title, inventory.inventory_id FROM sakila.film
RIGHT JOIN sakila.inventory ON film.film_id = inventory.film_id
WHERE film.film_id BETWEEN 13 AND 15; 

SELECT film.film_id, film.title, inventory.inventory_id FROM sakila.film
LEFT JOIN sakila.inventory ON film.film_id = inventory.film_id
WHERE film.film_id BETWEEN 13 AND 15; 



-- VIEW

CREATE VIEW customer_payment AS
SELECT customer_id, SUM(amount) AS total_payment
FROM sakila.payment
GROUP BY customer_id;

SELECT customer_id, total_payment
FROM customer_payment
WHERE total_payment > 100;



-- TRANSACTIONS

-- Lets create a table first

USE employee;

CREATE TABLE bank_account (
	account_id INT PRIMARY KEY,
    account_name VARCHAR(50),
    balance DECIMAL(10, 2)
);


INSERT INTO bank_account (account_id, account_name, balance) VALUES
("101", "Alice", 1000),
("102", "Bob", 1500),
("103", "Charlie", 2000),
("104", "David", "2500");


-- Let's transfer $300 from Alice's to Bob's acccount

START TRANSACTION;
-- Deduct $300 from Alice's account
UPDATE bank_account
SET balance = balance -300
WHERE account_id = 101;
-- Add $300 to Bob's account
UPDATE bank_account
SET balance=balance+300
WHERE account_id = 102;

COMMIT;


-- Savepoints can save time sometimes
UPDATE bank_account SET balance=balance-100 WHERE account_id=104;
SAVEPOINT meowmeow;

UPDATE bank_account SET balance=balance+100 WHERE account_id=101;

RELEASE SAVEPOINT meowmeow;
COMMIT;


-- INDEX

CREATE INDEX idx_name_balance ON bank_account(account_name, balance);
SELECT * FROM bank_account WHERE account_name="Charlie" AND balance >=1500;

DROP INDEX idx_name_balance ON bank_account;