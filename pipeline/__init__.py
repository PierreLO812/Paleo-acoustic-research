"""Paleo multi-agent pipeline: HYPOTHESIS → BUILD_AND_SIM → PEER_REVIEW → RETRY_LOOP."""

__all__ = ["run_pipeline"]


def __getattr__(name: str):
    if name == "run_pipeline":
        from pipeline.orchestrator import run_pipeline

        return run_pipeline
    raise AttributeError(name)
