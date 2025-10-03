import typing as t

from casp.services.api import schemas
from casp.services.api.services.annotation_engines.model_interpretation_engine.strategies.base import (
    ModelInterpretationStrategyInterface,
)
from casp.services.api.services.annotation_engines.shared_resources import CellOntologyResource
from casp.services.model_inference.schemas import ClassificationModelOutput


class OntologyAncestorUpliftStrategy(ModelInterpretationStrategyInterface):

    def __init__(
        self,
        prune_threshold: float,
        cell_ontology_resource: t.Optional[CellOntologyResource] = None,
    ):
        self.cell_ontology_resource = cell_ontology_resource or CellOntologyResource()
        self.prune_threshold = prune_threshold

    def _process_ancestor_uplifting(
        self,
        query_cell_id: str,
        probabilities,
        cl_classes_from_model: t.Sequence[str],
        cl_classes_from_model_set: t.Set[str],
    ):
        scores_dict = {k: 0 for k in self.cell_ontology_resource.ancestors_dictionary.keys()}

        # CAS needs those values for ontology aware response, this way it will be fully compatible.
        total_weight = 0
        total_neighbors = 0
        total_neighbors_unrecognized = 0

        for i, cl_class in enumerate(cl_classes_from_model):
            current_node_probability = probabilities[i]
            scores_dict[cl_class] = current_node_probability

            for ancestor in self.cell_ontology_resource.ancestors_dictionary[cl_class]:
                if ancestor in cl_classes_from_model_set:
                    break

                scores_dict[ancestor] = max(current_node_probability, scores_dict[ancestor])

        if self.prune_threshold > 0.0:
            scores_dict = {k: v for k, v in scores_dict.items() if v >= self.prune_threshold}

        # Normalize values
        max_score = max(scores_dict.values())
        scores_dict = {k: v / max_score for k, v in scores_dict.items()}

        annotation_summary = [
            schemas.AnnotationSummaryOntologyAware(
                score=score,
                cell_type_ontology_term_id=cell_type_ontology_term_id,
                cell_type=self.cell_ontology_resource.ontology_term_id_to_name_dict[cell_type_ontology_term_id],
            )
            for cell_type_ontology_term_id, score in scores_dict.items()
        ]
        return schemas.QueryCellAnnotationOntologyAware(
            query_cell_id=query_cell_id,
            matches=annotation_summary,
            total_weight=total_weight,
            total_neighbors=total_neighbors,
            total_neighbors_unrecognized=total_neighbors_unrecognized,
        )

    def summarize(self, model_output: ClassificationModelOutput) -> t.List[schemas.QueryCellAnnotationAbstract]:
        cl_classes_from_model = [x.replace(":", "_") for x in model_output.labels]
        cl_classes_from_model_set = set(cl_classes_from_model)  # precompute once to avoid repeating
        result = []

        for query_cell_id, probabilities in zip(model_output.sample_ids, model_output.probabilities):
            cell_annotation = self._process_ancestor_uplifting(
                query_cell_id=query_cell_id,
                probabilities=probabilities,
                cl_classes_from_model=cl_classes_from_model,
                cl_classes_from_model_set=cl_classes_from_model_set,
            )
            result.append(cell_annotation)

        return result
