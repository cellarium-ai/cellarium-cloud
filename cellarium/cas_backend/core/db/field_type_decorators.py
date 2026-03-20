import typing as t

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB, UUID


class SQLiteCompatibleJSONBType(sa.types.TypeDecorator):
    """
    Custom JSONB type to replace JSONB with Text or JSON when running tests on SQLite.
    """

    impl = JSONB

    def load_dialect_impl(self, dialect: sa.engine.interfaces.Dialect):
        if dialect.name == "sqlite":
            return dialect.type_descriptor(sa.types.JSON())  # Use JSON type for SQLite
        return dialect.type_descriptor(JSONB)

    def process_bind_param(self, value: t.Any, dialect: sa.engine.interfaces.Dialect):
        return value

    def process_result_value(self, value: t.Any, dialect: sa.engine.interfaces.Dialect):
        return value


class SQLiteCompatibleUUIDType(sa.types.TypeDecorator):
    """
    Custom UUID type for compatibility with tests on SQLite.
    """

    impl = UUID

    def load_dialect_impl(self, dialect: sa.engine.interfaces.Dialect):
        if dialect.name == "sqlite":
            return dialect.type_descriptor(sa.types.String(36))  # Store UUID as string in SQLite
        return dialect.type_descriptor(UUID)

    def process_bind_param(self, value: t.Any, dialect: sa.engine.interfaces.Dialect):
        if value is None:
            return None
        if dialect.name == "sqlite":
            return str(value)  # Convert UUID to string for SQLite
        return value

    def process_result_value(self, value: t.Any, dialect: sa.engine.interfaces.Dialect):
        if value is None:
            return None
        if dialect.name == "sqlite":
            return UUID(value)  # Convert string back to UUID
        return value
