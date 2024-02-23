import pickle
import typing as t

import redis

from casp.services import settings


class RedisSerializer:
    @staticmethod
    def dumps(obj: t.Any) -> bytes | int:
        """
        Serialize object to bytes.

        :param obj: Object to serialize. Can be anything that can be pickled.

        :return: Serialized object as bytes.
        """
        if type(obj) is int:
            # import django.core.cache.backends.base
            return obj
        return pickle.dumps(obj, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def loads(data: t.Any) -> t.Any:
        """
        Deserialize bytes to object.

        :param data: Serialized object as bytes.

        :return: Deserialized object.
        """
        # print(data)
        try:
            return int(data)
        except ValueError:
            return pickle.loads(data)


class RedisCache:
    """
    Redis cache client.
    """

    def __init__(self, cache_key_prefix: str = "main"):
        """
        Initialize Redis cache client.

        :param cache_key_prefix: Prefix to use for all cache keys.
        """
        self.serializer = RedisSerializer()
        self.cache_key_prefix = cache_key_prefix

    def _compose_key(self, key: str) -> str:
        """
        Compose cache key with prefix.

        :param key: Name of key to compose with prefix.

        :return: Composed key with prefix.
        """
        return f"{self.cache_key_prefix}:{key}"

    @staticmethod
    def get_client() -> redis.Redis:
        """
        Get Redis client.

        :return: Redis client.
        """
        return redis.Redis(
            host=settings.REDIS_HOST, port=settings.REDIS_PORT, db=settings.REDIS_DB, password=settings.REDIS_PASSWORD
        )

    def get(self, key: str, default: t.Optional[t.Any] = None) -> t.Any:
        """
        Get value from cache by key. Deserialize value if it exists. Return default if value does not exist.

        :param key: Name of key to get value from cache.
        :param default: Default value to return if key does not exist in cache.

        :return: Deserialized value from cache if it exists, default value otherwise.
        """
        client = self.get_client()
        _redis_key = self._compose_key(key)
        value = client.get(name=_redis_key)
        return default if value is None else self.serializer.loads(value)

    def set(self, key: str, value: t.Any, timeout: t.Optional[int] = settings.REDIS_CACHE_DEFAULT_TTL) -> None:
        """
        Set value in cache by key. Serialize value before setting in cache. Set default timeout if ttl is not provided.
        :param key: Name of key to set value in cache.
        :param value: Value to set in cache.
        :param timeout: Time to live for value in cache. If not provided, default timeout is used. If ``None``, no
            timeout is set.
        """
        client = self.get_client()
        _redis_key = self._compose_key(key)
        value = self.serializer.dumps(value)
        client.set(name=_redis_key, value=value, ex=timeout)

    def delete(self, key: str) -> bool:
        """
        Delete value from cache by key.

        :param key: Name of key to delete value from cache.

        :return: True if key exists in cache and value is deleted, False otherwise.
        """
        client = self.get_client()
        _redis_key = self._compose_key(key)
        return client.delete(_redis_key) > 0

    def exists(self, key: str) -> bool:
        """
        Check if key exists in cache.

        :param key: Name of key to check if it exists in cache.

        :return: True if key exists in cache, False otherwise.
        """
        client = self.get_client()
        _redis_key = self._compose_key(key)
        return client.exists(_redis_key) > 0

    def get_many(self, keys: t.List[str]) -> t.Dict[str, t.Any]:
        """
        Get multiple values from cache by keys. Deserialize values if they exist.

        :param keys: List of keys to get values from cache.

        :return: Dictionary of deserialized values from cache.
        """
        client = self.get_client()
        _redis_keys = [self._compose_key(key) for key in keys]
        redis_response = client.mget(_redis_keys)
        return {k: self.serializer.loads(v) for k, v in zip(keys, redis_response) if v is not None}

    def set_many(self, data: t.Dict[str, t.Any], timeout: t.Optional[int] = settings.REDIS_CACHE_DEFAULT_TTL) -> None:
        """
        Set multiple values in cache. Serialize values before setting in cache.

        :param data: Dictionary of key-value pairs to set in cache.
        :param timeout: Time to live for values in cache. If not provided, default timeout is used. If ``None``, no
            timeout is set. It is recommended to avoid using ``None`` as it can lead to memory leaks.
        """
        client = self.get_client()
        _redis_data = {self._compose_key(k): self.serializer.dumps(v) for k, v in data.items()}
        pipeline = client.pipeline()
        pipeline.mset(_redis_data)

        if timeout is not None:
            # Setting timeout for each key as redis does not support timeout
            # with mset().
            for key, _ in _redis_data.items():
                pipeline.expire(key, timeout)

        pipeline.execute()

    def exists_many(self, keys: t.List[str]) -> t.List[bool]:
        """
        Check if keys exist in cache.

        :param keys: List of keys to check if they exist in cache.

        :return: List of booleans indicating if keys exist in cache.
        """
        client = self.get_client()
        _redis_keys = [self._compose_key(key) for key in keys]
        return [client.exists(key) > 0 for key in _redis_keys]
