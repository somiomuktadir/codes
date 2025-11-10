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