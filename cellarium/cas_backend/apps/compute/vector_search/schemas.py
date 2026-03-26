from pydantic import BaseModel, Field


class MatchResult(BaseModel):
    """
    MatchResult represents the result of a match request to the matching service.
    """

    class Neighbor(BaseModel):
        """
        Neighbor represents a neighbor of a feature vector (e.g. a cell) in the matching service.
        """

        cas_cell_index: str = Field(default=None)
        distance: float = Field(default=None)

    class NearestNeighbors(BaseModel):
        """
        NearestNeighbors represents the result of a nearest neighbors request to the matching service. There should be 1
        per query that was sent.
        """

        neighbors: list["MatchResult.Neighbor"] = Field(default=[])

    matches: list["MatchResult.NearestNeighbors"] = Field(default=[])

    def concat(self, other: "MatchResult") -> "MatchResult":
        """
        Concatenate the results of two match requests.

        :param other: The other match result to concatenate.

        :return: The concatenated match result as a new instance of MatchResult.
        """
        return MatchResult(matches=self.matches + other.matches)

    def get_unique_ids(self) -> list[int]:
        """
        Get unique cas cell ids across all neighbors in querying cells in :class:`MatchResult`

        :return: A list with unique indexes of neighbors
        """
        return list(
            {int(neighbor.cas_cell_index) for query_neighbors in self.matches for neighbor in query_neighbors.neighbors}
        )


MatchResult.NearestNeighbors.update_forward_refs()
MatchResult.Neighbor.update_forward_refs()
