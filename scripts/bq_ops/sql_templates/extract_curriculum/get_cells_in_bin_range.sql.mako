<%
    from casp.data_manager.sql import mako_helpers as mh

    select_columns_processed = mh.add_cell_info_required_columns(select_columns)
    # Remove column aliases because after creating temporary tables all the columns are
    # part of the single table `c`
    select_columns_processed = mh.remove_leading_alias(select_columns_processed)
%>
${mh.select(select_columns_processed)}
from `${project}.${dataset}.${extract_table_prefix}__extract_cell_info` c
where extract_bin between ${start_bin} and ${end_bin}
order by farm_finger