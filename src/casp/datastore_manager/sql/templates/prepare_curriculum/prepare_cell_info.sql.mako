<%
    from casp.datastore_manager.sql import mako_helpers as mh

    select_column_names_processed = mh.add_cell_info_required_columns(select_columns)
    # Remove column aliases because after creating temporary tables all the columns are
    # part of the single table `c`
    select_column_names_processed = mh.remove_leading_alias(select_column_names_processed)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_cell_info`
partition by range_bucket(extract_bin, generate_array(0, ${partition_bin_count}, ${partition_size}))
cluster by extract_bin
as
${mh.select(select_column_names_processed)},
    cast(floor((row_number() over () - 1) / ${extract_bin_size}) as int) as extract_bin
from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info_randomized` c