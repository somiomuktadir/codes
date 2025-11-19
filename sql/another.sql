/*
 * SQL: Indexing and Transactions
 * Practical exercises covering essential indexing and transaction operations
 */

-- ============================================
-- PART 1: CREATE TABLES
-- ============================================

CREATE SCHEMA IF NOT EXISTS banking_demo;
USE banking_demo;

-- Accounts table
CREATE TABLE accounts (
    account_id INT PRIMARY KEY AUTO_INCREMENT,
    account_holder VARCHAR(100) NOT NULL,
    account_type ENUM('Checking', 'Savings', 'Business') DEFAULT 'Checking',
    balance DECIMAL(15, 2) NOT NULL DEFAULT 0.00,
    branch_code VARCHAR(10),
    city VARCHAR(50),
    created_date DATE DEFAULT (CURRENT_DATE)
);

-- Transaction log table
CREATE TABLE transaction_log (
    transaction_id BIGINT PRIMARY KEY AUTO_INCREMENT,
    from_account INT,
    to_account INT,
    transaction_type ENUM('Deposit', 'Withdrawal', 'Transfer') NOT NULL,
    amount DECIMAL(15, 2) NOT NULL,
    transaction_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    description VARCHAR(255),
    FOREIGN KEY (from_account) REFERENCES accounts(account_id),
    FOREIGN KEY (to_account) REFERENCES accounts(account_id)
);

-- Insert sample data
INSERT INTO accounts (account_holder, account_type, balance, branch_code, city) VALUES
('Alice Johnson', 'Checking', 5000.00, 'BR001', 'New York'),
('Bob Smith', 'Savings', 15000.00, 'BR001', 'New York'),
('Charlie Brown', 'Business', 50000.00, 'BR002', 'Los Angeles'),
('Diana Prince', 'Checking', 8000.00, 'BR002', 'Los Angeles'),
('Edward Norton', 'Savings', 25000.00, 'BR003', 'Chicago');


-- ============================================
-- PART 2: INDEXING BASICS
-- ============================================

-- 2.1. Single Column Index
-- Speeds up queries filtering by account_holder
CREATE INDEX idx_account_holder ON accounts(account_holder);

SELECT * FROM accounts WHERE account_holder = 'Alice Johnson';

-- 2.2. Composite Index (Multi-Column)
-- Optimizes queries filtering by multiple columns
CREATE INDEX idx_city_branch ON accounts(city, branch_code);

SELECT * FROM accounts WHERE city = 'New York' AND branch_code = 'BR001';

-- 2.3. Unique Index
-- Ensures uniqueness and speeds up lookups
CREATE UNIQUE INDEX idx_holder_branch ON accounts(account_holder, branch_code);

-- 2.4. Index on Foreign Keys
-- Speeds up JOIN operations
CREATE INDEX idx_from_account ON transaction_log(from_account);
CREATE INDEX idx_to_account ON transaction_log(to_account);

-- 2.5. View Indexes on a Table
SHOW INDEXES FROM accounts;

-- 2.6. Check if Index is Used (EXPLAIN)
EXPLAIN SELECT * FROM accounts WHERE city = 'New York';

-- 2.7. Drop an Index
-- DROP INDEX idx_account_holder ON accounts;


-- ============================================
-- PART 3: BASIC TRANSACTIONS
-- ============================================

-- 3.1. Simple Deposit
START TRANSACTION;

UPDATE accounts 
SET balance = balance + 1000.00 
WHERE account_id = 1;

INSERT INTO transaction_log (to_account, transaction_type, amount, description)
VALUES (1, 'Deposit', 1000.00, 'Cash deposit');

COMMIT;


-- 3.2. Withdrawal with Validation
START TRANSACTION;

UPDATE accounts 
SET balance = balance - 500.00 
WHERE account_id = 2 AND balance >= 500.00;

INSERT INTO transaction_log (from_account, transaction_type, amount, description)
VALUES (2, 'Withdrawal', 500.00, 'ATM withdrawal');

COMMIT;


-- 3.3. Transaction Rollback
START TRANSACTION;

UPDATE accounts 
SET balance = balance - 10000.00 
WHERE account_id = 3;

-- Oops, wrong amount! Undo the transaction
ROLLBACK;


-- ============================================
-- PART 4: MONEY TRANSFER
-- ============================================

-- 4.1. Transfer Between Accounts
START TRANSACTION;

-- Deduct from sender (account 1)
UPDATE accounts 
SET balance = balance - 250.00 
WHERE account_id = 1 AND balance >= 250.00;

-- Add to receiver (account 2)
UPDATE accounts 
SET balance = balance + 250.00 
WHERE account_id = 2;

-- Log the transfer
INSERT INTO transaction_log (from_account, to_account, transaction_type, amount, description)
VALUES (1, 2, 'Transfer', 250.00, 'Transfer between accounts');

COMMIT;

-- Verify the transfer
SELECT account_id, account_holder, balance 
FROM accounts 
WHERE account_id IN (1, 2);


-- ============================================
-- PART 5: SAVEPOINTS
-- ============================================

-- 5.1. Using Savepoints for Partial Rollback
START TRANSACTION;

-- Operation 1: Deposit
UPDATE accounts SET balance = balance + 500 WHERE account_id = 3;
SAVEPOINT after_deposit;

-- Operation 2: Apply fee
UPDATE accounts SET balance = balance - 10 WHERE account_id = 3;
SAVEPOINT after_fee;

-- Operation 3: Withdrawal
UPDATE accounts SET balance = balance - 200 WHERE account_id = 3;

-- Undo only the withdrawal, keep deposit and fee
ROLLBACK TO SAVEPOINT after_fee;

COMMIT;


-- ============================================
-- PART 6: MULTIPLE UPDATES IN ONE TRANSACTION
-- ============================================

-- 6.1. Batch Update - Apply Interest to All Savings Accounts
START TRANSACTION;

-- Add 0.5% interest to all savings accounts
UPDATE accounts 
SET balance = balance * 1.005
WHERE account_type = 'Savings';

-- Log interest for each account
INSERT INTO transaction_log (to_account, transaction_type, amount, description)
SELECT 
    account_id,
    'Deposit',
    balance * 0.005,
    'Monthly interest'
FROM accounts
WHERE account_type = 'Savings';

COMMIT;


-- ============================================
-- CLEANUP (Optional)
-- ============================================

-- Drop tables if needed
-- DROP TABLE IF EXISTS transaction_log;
-- DROP TABLE IF EXISTS accounts;
-- DROP SCHEMA IF EXISTS banking_demo;


-- ============================================
-- KEY TAKEAWAYS
-- ============================================

/*
INDEXING:
- Create indexes on columns used in WHERE, JOIN, ORDER BY
- Use composite indexes for multi-column queries
- Check index usage with EXPLAIN
- Don't over-index (slows INSERT/UPDATE)

TRANSACTIONS:
- Always use START TRANSACTION and COMMIT
- Use ROLLBACK to undo changes
- Use savepoints for partial rollbacks
- Keep transactions short and focused
- Validate data before committing
*/
