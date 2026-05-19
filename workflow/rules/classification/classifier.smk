# Classification rule selector
# Chooses a single classifier backend per run.
# Uses shared helpers from rules/common.smk.

CLASSIFICATION_CONFIG = get_module_config(config, "classification")

engine = get_classification_engine(config)

if engine == "rdp":
    include: "rdp_classifier.smk"
elif engine == "sintax":
    include: "sintax_classifier.smk"
else:
    raise ValueError(
        f"Unknown classification engine={engine!r}. Supported engines: rdp, sintax"
    )
