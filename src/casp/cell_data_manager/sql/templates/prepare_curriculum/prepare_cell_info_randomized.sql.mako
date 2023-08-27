<%
    from casp.cell_data_manager.sql import mako_helpers

    select_columns_processed = mako_helpers.add_cell_info_required_columns(select_columns)
    select_columns_parsed = mako_helpers.parse_column_names(select_columns_processed)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_cell_info_randomized`
as
select ${select_columns_parsed}
from `${project}.${dataset}.cas_cell_info` c
join `${project}.${dataset}.cas_ingest_info` i on (i.cas_ingest_id = c.cas_ingest_id)
where ${mako_helpers.parse_where_body(filter_statements)}
order by farm_fingerprint(cast(cas_cell_index + ${random_seed_offset} as string))