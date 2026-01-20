# Classification rule selector
# Chooses a single classifier backend per run.
# Prefer `config["modules"]["classification_engine"]` over legacy `config["classification"]["engine"]`.

def get_module_config(config, module_name):
    """Safely get module configuration, handling boolean values."""
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        module_config = modules.get(module_name, {})
        if isinstance(module_config, dict):
            return module_config
    fallback = config.get(module_name, {})
    return fallback if isinstance(fallback, dict) else {}


CLASSIFICATION_CONFIG = get_module_config(config, "classification")

def get_classification_engine(config, classification_config):
    modules = config.get("modules", {})
    if isinstance(modules, dict):
        value = modules.get("classification_engine")
        if value is None:
            value = modules.get("classifier")
        nested = modules.get("classification")
        if value is None and isinstance(nested, dict):
            value = nested.get("classification_engine") or nested.get("classifier") or nested.get("engine")
        if value is not None:
            return str(value).strip().lower()

    value = classification_config.get("engine")
    if value is not None:
        return str(value).strip().lower()
    return "rdp"


engine = get_classification_engine(config, CLASSIFICATION_CONFIG)

if engine == "rdp":
    include: "rdp_classifier.smk"
elif engine == "sintax":
    include: "sintax_classifier.smk"
else:
    raise ValueError(
        f"Unknown classification engine={engine!r}. Supported engines: rdp, sintax"
    )
