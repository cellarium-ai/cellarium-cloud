<%
    from casp.cell_data_manager.sql import mako_helpers

    select_columns_processed = mako_helpers.add_cell_info_required_columns(select_columns)
    # Remove column aliases because after creating temporary tables all the columns are
    # part of the single table `c`
    select_columns_processed = mako_helpers.remove_leading_alias(select_columns_processed)
    select_columns_parsed = mako_helpers.parse_column_names(select_columns_processed)
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_cell_info`
partition by range_bucket(extract_bin, generate_array(0,${partition_bin_count},${partition_size}))
cluster by extract_bin
as
select ${select_columns_parsed},
        cast(floor((row_number() over () - 1) / ${extract_bin_size}) as int) as extract_bin
from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info_randomized` c