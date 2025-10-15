import typing as t
from collections import defaultdict

import numpy as np
from casp.services.api import schemas
from casp.services.api.services.annotation_engines.model_interpretation_engine.strategies.base import (
    ModelInterpretationStrategyInterface,
)
from casp.services.api.services.annotation_engines.shared_resources import CellOntologyResource
from casp.services.model_inference.schemas import ClassificationModelOutput

ROOT_NODE_CL_LABEL = "CL_0000000"


class OntologyAncestorUpliftStrategy(ModelInterpretationStrategyInterface):
    def __init__(
        self,
        prune_threshold: float = 0.0,
        decay: float = 1.0,  # 1.0 = no decay; <1.0 attenuates farther ancestors
        stop_at_model_ancestor: bool = True,  # True = stop climbing once you hit a model-covered ancestor
        top_k: t.Optional[int] = None,
        cell_ontology_resource: t.Optional[CellOntologyResource] = None,
    ):
        self.cell_ontology_resource = cell_ontology_resource or CellOntologyResource()
        self.prune_threshold = prune_threshold
        self.decay = decay
        self.stop_at_model_ancestor = stop_at_model_ancestor
        self.top_k = top_k

    def _process_ancestor_uplifting(
        self,
        query_cell_id: str,
        probabilities: np.ndarray,  # 1D array, aligned with model label order
        cl_classes_from_model: t.Sequence[str],  # "CL_XXXX"
        cl_classes_from_model_set: t.Set[str],
    ) -> schemas.QueryCellAnnotationOntologyAware:

        scores = defaultdict(float)

        for i, cl_term in enumerate(cl_classes_from_model):
            p = float(probabilities[i])
            # set the class itself
            scores[cl_term] = max(scores[cl_term], p)
            # scores[cl_term] = p

            # uplift ancestors
            for depth, anc in enumerate(self.cell_ontology_resource.ancestors_dictionary.get(cl_term, []), start=1):
                if self.stop_at_model_ancestor and anc in cl_classes_from_model_set:
                    break  # we already updated it; don't go higher if you choose this policy
                w = p * (self.decay**depth)
                # scores[anc] = max(scores[anc], w)
                scores[anc] = scores[anc] + w

        # prune
        if self.prune_threshold > 0.0:
            scores = {k: v for k, v in scores.items() if v >= self.prune_threshold}

        # ensure non-empty (keep top1 if everything pruned)
        if not scores:
            j = int(np.argmax(probabilities))
            scores = {cl_classes_from_model[j]: float(probabilities[j])}

        # sort + optional top-k
        items = sorted(scores.items(), key=lambda kv: kv[1], reverse=True)
        if self.top_k:
            items = items[: self.top_k]

        # build payload; be tolerant if a name is missing
        matches = [
            schemas.AnnotationSummaryOntologyAware(
                score=score,
                cell_type_ontology_term_id=term_id,
                cell_type=self.cell_ontology_resource.ontology_term_id_to_name_dict.get(term_id, term_id),
            )
            for term_id, score in items
        ]

        # CAS compat fields; keep your placeholders
        return schemas.QueryCellAnnotationOntologyAware(
            query_cell_id=query_cell_id,
            matches=matches,
            total_weight=0,
            total_neighbors=0,
            total_neighbors_unrecognized=0,
        )

        # root_node_cell_type = self.cell_ontology_resource.ontology_term_id_to_name_dict.get(
        #     ROOT_NODE_CL_LABEL, ROOT_NODE_CL_LABEL
        # )
        # root_node_annotation = schemas.AnnotationSummaryOntologyAware(
        #     score=1, cell_type_ontology_term_id=ROOT_NODE_CL_LABEL, cell_type=root_node_cell_type
        # )
        # matches = [root_node_annotation]
        # for cl_class, probability in zip(cl_classes_from_model, probabilities):
        #     if probability >= self.prune_threshold:
        #         matches.append(
        #             schemas.AnnotationSummaryOntologyAware(
        #                 score=probability,
        #                 cell_type_ontology_term_id=cl_class,
        #                 cell_type=self.cell_ontology_resource.ontology_term_id_to_name_dict.get(cl_class, cl_class),
        #             )
        #         )

        return schemas.QueryCellAnnotationOntologyAware(
            query_cell_id=query_cell_id,
            matches=matches,
            total_weight=0,
            total_neighbors=0,
            total_neighbors_unrecognized=0,
        )

    def summarize(self, model_output: ClassificationModelOutput) -> t.List[schemas.QueryCellAnnotationAbstract]:
        # normalize your IDs once; make sure your resource uses the same convention
        model_terms = [x.replace(":", "_") for x in model_output.labels]
        model_set = set(model_terms)
        out = []
        for query_cell_id, probs in zip(model_output.sample_ids, model_output.probabilities):
            out.append(
                self._process_ancestor_uplifting(
                    query_cell_id=query_cell_id,
                    probabilities=probs,
                    cl_classes_from_model=model_terms,
                    cl_classes_from_model_set=model_set,
                )
            )
        return out
