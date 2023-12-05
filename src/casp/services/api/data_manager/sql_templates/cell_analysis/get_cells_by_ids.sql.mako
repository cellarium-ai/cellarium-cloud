<%
    from casp.data_manager.sql import mako_helpers as mh
%>
${mh.select(select_columns)}
from `${project}.${dataset}.cas_cell_info` ci
${mh.where(filter_statements)}