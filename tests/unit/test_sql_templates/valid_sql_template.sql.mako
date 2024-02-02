<%!
    from casp.data_manager.sql import mako_helpers as mh
%>
${mh.select(select_columns)}
from `${project}.${dataset}.test_table_name` ttn
join `${project}.${dataset}.test_join_table` tjt on (ttn.id = tjt.id)
${mh.where(filter_statements)}
group by 1
limit 10