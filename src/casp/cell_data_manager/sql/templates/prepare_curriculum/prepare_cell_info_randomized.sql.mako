<%
    from casp.cell_data_manager.sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_cell_info_randomized`
as
${mh.select(select_column_names_processed)}
from `${project}.${dataset}.cas_cell_info` c
join `${project}.${dataset}.cas_ingest_info` i on (i.cas_ingest_id = c.cas_ingest_id)
${mh.where(filter_statements)}
order by farm_fingerprint(cast(cas_cell_index + ${random_seed_offset} as string))