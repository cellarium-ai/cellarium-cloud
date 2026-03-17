update `${project}.${dataset}.cas_cell_info` ci
set bq_row_number = t.row_num
from (
    select
        row_number() over (order by cas_cell_index) as row_num,
        cas_cell_index
    from `${project}.${dataset}.cas_cell_info`
) as t
where ci.cas_cell_index = t.cas_cell_index
