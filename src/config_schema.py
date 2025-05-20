# src/config_schema.py

config_schema = {
	"path": {"type": "string", "required": True},
	"trajectory_file": {"type": "string", "required": True},
	"topology_file": {"type": "string", "required": True},
	"number_of_molecules": {"type": "integer", "min": 2, "required": True},
	"center_atom": {"type": "string", "min": 1, "required": True},
	"micelle_residue_atom": {"type": "list", "required": True},
	"ion_radius": {"type": "integer", "min": 2, "required": True},
	"cutoff_distances": {"type": "list","schema":{"type": "float"}, "required": True},
	"number_cpu": {"type": "integer", "min": 1, "required": False,"default":1},
	"name":{"type": "string", "required": True}

}
