-- --- Exercise 1: Creating Tables and Adding Constraints ---

-- 1. Create the employee table.
CREATE TABLE employee (
    employee_id SMALLINT UNSIGNED PRIMARY KEY,
    first_name VARCHAR(20) NOT NULL,
    last_name VARCHAR(20) NOT NULL,
    job_title VARCHAR(20),
    salary DECIMAL(10,2),
    department VARCHAR(20)
);

-- 2. Create the project table, with a foreign key constraint.
CREATE TABLE project (
    project_id SMALLINT UNSIGNED PRIMARY KEY,
    project_name VARCHAR(50),
    employee_id SMALLINT UNSIGNED,
    FOREIGN KEY (employee_id) REFERENCES employee(employee_id)
);

-- --- Exercise 2: Modifying and Inserting Data ---

-- 1. Make employee_id in the employee table an AUTO_INCREMENT column.
ALTER TABLE employee MODIFY employee_id SMALLINT UNSIGNED AUTO_INCREMENT PRIMARY KEY;

-- 2. Insert the sample data into the employee table.
INSERT INTO employee (employee_id, first_name, last_name, job_title, salary, department) VALUES
(NULL, 'Alice', 'Smith', 'Developer', 60000, 'IT'),
(NULL, 'Bob', 'Brown', 'Manager', 75000, 'Sales');

-- --- Exercise 3: Querying Data ---

-- 1. Retrieve the first_name, last_name, and job_title of each employee.
SELECT first_name, last_name, job_title FROM employee;

-- 2. Retrieve employees with a salary greater than 65000.
SELECT * FROM employee WHERE salary > 65000;

-- 3. Retrieve the project_name and employee_id from the project table for a specific employee_id (e.g., employee_id = 1).
SELECT project_name, employee_id FROM project WHERE employee_id = 1;

-- --- Exercise 4: Updating Data ---

-- Update the salary of the employee with employee_id = 1 to 62000.
UPDATE employee SET salary = 62000 WHERE employee_id = 1;

-- --- Exercise 5: Deleting Data ---

-- Delete the project record where project_id = 1.
DELETE FROM project WHERE project_id = 1;

-- --- Exercise 6: Dropping Tables ---

-- Drop the project table.
DROP TABLE project;

-- Drop the employee table.
DROP TABLE employee;