-- 01_advanced_sakila_queries.sql
-- Advanced queries for the Sakila sample database.

-- 1. List the top 5 customers who have rented the most movies, 
-- along with their total spending.
SELECT 
    c.customer_id,
    c.first_name,
    c.last_name,
    COUNT(r.rental_id) AS total_rentals,
    SUM(p.amount) AS total_spent
FROM customer c
JOIN rental r ON c.customer_id = r.customer_id
JOIN payment p ON r.rental_id = p.rental_id
GROUP BY c.customer_id
ORDER BY total_rentals DESC
LIMIT 5;

-- 2. Find the film categories that have an average rental duration 
-- greater than the overall average rental duration of all films.
SELECT 
    c.name AS category_name,
    AVG(f.rental_duration) AS avg_duration
FROM category c
JOIN film_category fc ON c.category_id = fc.category_id
JOIN film f ON fc.film_id = f.film_id
GROUP BY c.name
HAVING AVG(f.rental_duration) > (
    SELECT AVG(rental_duration) FROM film
);

-- 3. List all actors who have appeared in a movie with 'ACADEMY' in the title.
SELECT DISTINCT
    a.first_name,
    a.last_name
FROM actor a
JOIN film_actor fa ON a.actor_id = fa.actor_id
JOIN film f ON fa.film_id = f.film_id
WHERE f.title LIKE '%ACADEMY%';

-- 4. Revenue per store per month for the year 2005.
SELECT 
    s.store_id,
    EXTRACT(MONTH FROM p.payment_date) AS payment_month,
    SUM(p.amount) AS monthly_revenue
FROM store s
JOIN staff st ON s.store_id = st.store_id
JOIN payment p ON st.staff_id = p.staff_id
WHERE EXTRACT(YEAR FROM p.payment_date) = 2005
GROUP BY s.store_id, payment_month
ORDER BY s.store_id, payment_month;

-- 5. Find customers who have rented movies from both stores.
SELECT customer_id, first_name, last_name
FROM customer
WHERE customer_id IN (
    SELECT DISTINCT r.customer_id
    FROM rental r
    JOIN inventory i ON r.inventory_id = i.inventory_id
    WHERE i.store_id = 1
)
AND customer_id IN (
    SELECT DISTINCT r.customer_id
    FROM rental r
    JOIN inventory i ON r.inventory_id = i.inventory_id
    WHERE i.store_id = 2
);
