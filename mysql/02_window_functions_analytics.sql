-- 02_window_functions_analytics.sql
-- Analytics queries using Window Functions and CTEs.

-- 1. Rank customers by total spending using DENSE_RANK.
WITH CustomerSpending AS (
    SELECT 
        c.customer_id,
        c.first_name,
        c.last_name,
        SUM(p.amount) AS total_amount
    FROM customer c
    JOIN payment p ON c.customer_id = p.customer_id
    GROUP BY c.customer_id
)
SELECT 
    RANK() OVER (ORDER BY total_amount DESC) AS rank_val,
    DENSE_RANK() OVER (ORDER BY total_amount DESC) AS dense_rank_val,
    first_name,
    last_name,
    total_amount
FROM CustomerSpending
ORDER BY total_amount DESC
LIMIT 20;

-- 2. Calculate the running total of monthly revenue for each store.
WITH MonthlyRevenue AS (
    SELECT 
        s.store_id,
        DATE_FORMAT(p.payment_date, '%Y-%m') AS payment_month,
        SUM(p.amount) AS revenue
    FROM store s
    JOIN staff st ON s.store_id = st.store_id
    JOIN payment p ON st.staff_id = p.staff_id
    GROUP BY s.store_id, payment_month
)
SELECT 
    store_id,
    payment_month,
    revenue,
    SUM(revenue) OVER (PARTITION BY store_id ORDER BY payment_month) AS running_total
FROM MonthlyRevenue;

-- 3. Find the difference in rental count for each customer compared to the previous customer 
-- (when sorted by rental count).
WITH CustomerRentals AS (
    SELECT 
        customer_id,
        COUNT(rental_id) AS rental_count
    FROM rental
    GROUP BY customer_id
)
SELECT 
    customer_id,
    rental_count,
    LAG(rental_count, 1) OVER (ORDER BY rental_count DESC) AS prev_rental_count,
    rental_count - LAG(rental_count, 1) OVER (ORDER BY rental_count DESC) AS diff_from_prev
FROM CustomerRentals
ORDER BY rental_count DESC;

-- 4. Identify the top 3 films in each category by rental count.
WITH FilmRentals AS (
    SELECT 
        c.name AS category_name,
        f.title,
        COUNT(r.rental_id) AS rental_count,
        ROW_NUMBER() OVER (PARTITION BY c.name ORDER BY COUNT(r.rental_id) DESC) AS row_num
    FROM category c
    JOIN film_category fc ON c.category_id = fc.category_id
    JOIN film f ON fc.film_id = f.film_id
    JOIN inventory i ON f.film_id = i.film_id
    JOIN rental r ON i.inventory_id = r.inventory_id
    GROUP BY c.name, f.title
)
SELECT 
    category_name,
    title,
    rental_count
FROM FilmRentals
WHERE row_num <= 3;
