select
coalesce(sum(cells_processed), 0)
from (select distinct
        request_id,
        array_agg(event) events, -- aggregates all of the requests ids into an array
        max(cell_count) cells_processed -- presuming the start or end is 0 and the other is the cell count
    from users_useractivity ua
    where user_id = :id
        and finished_time >= :start_of_week
        and request_id is not null
    group by request_id
    ) requests
where not 'FAILED' = ANY(events);