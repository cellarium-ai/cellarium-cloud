<%
    from casp.cell_data_manager.sql import mako_helpers

    select_columns_processed = mako_helpers.add_cell_info_required_columns(select_columns)
    # Remove column aliases because after creating temporary tables all the columns are
    # part of the single table `c`
    select_columns_processed = mako_helpers.remove_leading_alias(select_columns_processed)
    select_columns_parsed = mako_helpers.parse_column_names(select_columns_processed)
%>
select ${select_columns_parsed}
from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info` c
where extract_bin between ${start_bin} and ${end_bin}