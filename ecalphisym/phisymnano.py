from typing import Any, Dict
from coffea.nanoevents.schemas import BaseSchema
import urllib.parse

def _build_record_array(
        name: str,
        name_mapping: Dict[str, str],
        contents: Dict[str, Any],
        record_name: str,
        class_name: str = None,
        size: int = None ) -> Dict[str, Any]:
    """Build a record array using the mapping we've got from the contents.

    Args:
        name_mapping (Dict[str, str]): The mapping of user variable to column in the contents
        contents (Dict[str, Any]): The contents of the array we are building into
    """
    items = {
        v_name: contents[col_name]["content"] if "content" in contents[col_name] else contents[col_name]
        for v_name, col_name in name_mapping.items()
    }
    record = {
        "class": "RecordArray",
        "contents": items,
        "form_key": urllib.parse.quote("!invalid," + name, safe=""),
        "parameters": {"__record__": record_name},
    }
    first = contents[next(iter(name_mapping.values()))]
    return {
        "class": class_name if class_name else first["class"],
        "size" : size,
        "offsets": None,
        "content": record,
        "form_key": first["form_key"],
        "parameters": {},
    }


class EcalPhiSymSchema(BaseSchema):
    """Build a schema using heuristics to imply a structure"""

    def __init__(self, base_form: Dict[str, Any]):
        """Create a EcalPhiSym schema.

        Notes:
            - Currently the repackaging of sumet_m1,sumet_m2,...,sumet_pN is done
              inside RecHit.sumet_v. This could probably be also done while generating
              the records here.

        Args:
            base_form (Dict[str, Any]): The base form of what we are going to generate a new schema (form) for.
        """
        super().__init__(base_form)

        # CMSSW to awkward mapping
        ecal_phisym_classes = {
            "EcalPhiSymEB" : { "record_name" : "RecHitEB" },
            "EcalPhiSymEE" : { "record_name" : "RecHitEE" },
            "EcalPhiSymInfo" : { "record_name" : "Info",
                                 "class_name" : "RegularArray",
                                 "size" : 1 }}
        
        # Get the collection names - anything with a common name before the "_".
        contents = self._form["contents"]
        collections = set(k.split("_")[0] for k in contents if "_" in k)
        
        output = {}
        for c_name in collections:
            mapping = {
                k.split(f"{c_name}_")[1].lower(): k for k in contents if k.startswith(f"{c_name}_")
            }
                        
            record = _build_record_array(
                c_name,
                mapping,
                contents,
                **ecal_phisym_classes[c_name]
            )

            record["parameters"].update({"collection_name": c_name})
            output[c_name] = record
            
        # Single items in the collection
        single_items = [k for k in contents if "_" not in k]
        for item_name in single_items:
            output[item_name] = contents[item_name]

        self._form["contents"] = output
        
    @property
    def behavior(self):
        """Behaviors necessary to implement this schema"""
        from coffea.nanoevents.methods import base
        from ecalphisym import rechit
        from ecalphisym import info

        behavior = {}
        behavior.update(base.behavior)
        behavior.update(rechit.behavior)
        behavior.update(info.behavior)
        return behavior
