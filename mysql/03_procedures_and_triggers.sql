-- 03_procedures_and_triggers.sql
-- Stored Procedures and Triggers for business logic and auditing.

DELIMITER //

-- 1. Stored Procedure: Get Customer Rental History
-- Returns the list of films rented by a specific customer.
CREATE PROCEDURE GetCustomerRentalHistory(IN cust_id INT)
BEGIN
    SELECT 
        f.title,
        r.rental_date,
        r.return_date
    FROM rental r
    JOIN inventory i ON r.inventory_id = i.inventory_id
    JOIN film f ON i.film_id = f.film_id
    WHERE r.customer_id = cust_id
    ORDER BY r.rental_date DESC;
END //

-- 2. Stored Procedure: Add New Film
-- Inserts a new film and links it to a language.
CREATE PROCEDURE AddNewFilm(
    IN p_title VARCHAR(255),
    IN p_description TEXT,
    IN p_release_year YEAR,
    IN p_language_id TINYINT,
    IN p_rental_duration TINYINT,
    IN p_rental_rate DECIMAL(4,2),
    IN p_length SMALLINT,
    IN p_replacement_cost DECIMAL(5,2),
    IN p_rating ENUM('G','PG','PG-13','R','NC-17')
)
BEGIN
    INSERT INTO film (
        title, description, release_year, language_id, 
        rental_duration, rental_rate, length, 
        replacement_cost, rating, last_update
    ) VALUES (
        p_title, p_description, p_release_year, p_language_id,
        p_rental_duration, p_rental_rate, p_length,
        p_replacement_cost, p_rating, NOW()
    );
END //

-- 3. Trigger: Audit Rental Updates
-- Logs changes to the rental table (e.g., when a movie is returned).
-- Assumes an 'audit_log' table exists.

/*
CREATE TABLE IF NOT EXISTS audit_log (
    log_id INT AUTO_INCREMENT PRIMARY KEY,
    table_name VARCHAR(50),
    action_type VARCHAR(10),
    record_id INT,
    log_timestamp DATETIME DEFAULT CURRENT_TIMESTAMP,
    details TEXT
);
*/

CREATE TRIGGER after_rental_update
AFTER UPDATE ON rental
FOR EACH ROW
BEGIN
    IF OLD.return_date IS NULL AND NEW.return_date IS NOT NULL THEN
        INSERT INTO audit_log (table_name, action_type, record_id, details)
        VALUES ('rental', 'RETURN', NEW.rental_id, CONCAT('Movie returned on ', NEW.return_date));
    END IF;
END //

-- 4. Function: Calculate Overdue Fees
-- Calculates potential overdue fees based on days late.
CREATE FUNCTION CalculateOverdueFees(p_rental_id INT) 
RETURNS DECIMAL(10,2)
DETERMINISTIC
BEGIN
    DECLARE v_days_overdue INT;
    DECLARE v_fee DECIMAL(10,2);
    DECLARE v_rental_rate DECIMAL(4,2);
    
    SELECT 
        DATEDIFF(NOW(), r.rental_date) - f.rental_duration,
        f.rental_rate
    INTO v_days_overdue, v_rental_rate
    FROM rental r
    JOIN inventory i ON r.inventory_id = i.inventory_id
    JOIN film f ON i.film_id = f.film_id
    WHERE r.rental_id = p_rental_id;
    
    IF v_days_overdue > 0 THEN
        SET v_fee = v_days_overdue * 1.00; -- $1.00 per day late fee
    ELSE
        SET v_fee = 0.00;
    END IF;
    
    RETURN v_fee;
END //

DELIMITER ;
