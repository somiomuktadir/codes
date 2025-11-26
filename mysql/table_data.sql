/*
 * SQL Table Data Management Exercises
 * Demonstrates DDL (Data Definition Language) and DML (Data Manipulation Language)
 * Covers: CREATE, ALTER, INSERT, SELECT, UPDATE, DELETE, DROP
 */

-- ============================================
-- EXERCISE 1: CREATING TABLES AND ADDING CONSTRAINTS
-- ============================================

-- 1. Create the employee table
CREATE TABLE employee (
    employee_id SMALLINT UNSIGNED PRIMARY KEY,
    first_name VARCHAR(20) NOT NULL,
    last_name VARCHAR(20) NOT NULL,
    job_title VARCHAR(20),
    salary DECIMAL(10, 2),
    department VARCHAR(20)
);

-- 2. Create the project table with a foreign key constraint
-- Foreign keys ensure referential integrity between tables
CREATE TABLE project (
    project_id SMALLINT UNSIGNED PRIMARY KEY,
    project_name VARCHAR(50),
    employee_id SMALLINT UNSIGNED,
    FOREIGN KEY (employee_id) REFERENCES employee(employee_id)
);


-- ============================================
-- EXERCISE 2: MODIFYING AND INSERTING DATA
-- ============================================

-- 1. Make employee_id an AUTO_INCREMENT column
-- AUTO_INCREMENT automatically generates unique IDs
ALTER TABLE employee 
MODIFY employee_id SMALLINT UNSIGNED AUTO_INCREMENT PRIMARY KEY;

-- 2. Insert sample data into the employee table
-- Using NULL for AUTO_INCREMENT columns lets MySQL generate the ID
INSERT INTO employee (first_name, last_name, job_title, salary, department) VALUES
('Alice', 'Smith', 'Developer', 60000, 'IT'),
('Bob', 'Brown', 'Manager', 75000, 'Sales');


-- ============================================
-- EXERCISE 3: QUERYING DATA
-- ============================================

-- 1. Retrieve specific columns for each employee
SELECT first_name, last_name, job_title 
FROM employee;

-- 2. Retrieve employees with a salary greater than 65000
SELECT * 
FROM employee 
WHERE salary > 65000;

-- 3. Retrieve projects for a specific employee
-- Replace employee_id = 1 with the desired employee ID
SELECT project_name, employee_id 
FROM project 
WHERE employee_id = 1;


-- ============================================
-- EXERCISE 4: UPDATING DATA
-- ============================================

-- Update the salary of a specific employee
-- Best practice: Always use WHERE clause to avoid updating all rows
UPDATE employee 
SET salary = 62000 
WHERE employee_id = 1;


-- ============================================
-- EXERCISE 5: DELETING DATA
-- ============================================

-- Delete a specific project record
-- Best practice: Always use WHERE clause to avoid deleting all rows
DELETE FROM project 
WHERE project_id = 1;


-- ============================================
-- EXERCISE 6: DROPPING TABLES
-- ============================================

-- Drop tables in correct order (child table first, then parent)
-- Must drop project table first due to foreign key constraint

-- Drop the project table
DROP TABLE IF EXISTS project;

-- Drop the employee table
DROP TABLE IF EXISTS employee;