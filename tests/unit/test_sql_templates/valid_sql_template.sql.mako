<%!
    from casp.cell_data_manager.sql import mako_helpers
%>
select ${mako_helpers.parse_column_names(select_columns)}
from `${project}.${dataset}.test_table_name` ttn
join `${project}.${dataset}.test_join_table` tjt on (ttn.id = tjt.id)
where ${mako_helpers.parse_where_body(filter_statements)}
group by 1
limit 10