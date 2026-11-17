"""
Orchestrates:
  HYPOTHESIS -> BUILD_AND_SIM -> PEER_REVIEW -> RETRY_LOOP (max 3)

Modes:
  - local: deterministic Python agents (default smoke path)
  - agents: emit prompts / hand off to Grokbot + Antigravity (T3)
"""

from __future__ import annotations

import argparse
import json
import uuid
from datetime import datetime, timezone

from pipeline.build_and_sim import run_build_and_sim
from pipeline.hypothesis import list_taxa, write_hypothesis_bundle
from pipeline.paths import OUTPUT
from pipeline.peer_review import normalize_legacy_state, review_run


def new_run_id(taxon_slug: str) -> str:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{taxon_slug}_{stamp}_{uuid.uuid4().hex[:6]}"


def run_pipeline(
    taxon_slug: str,
    *,
    max_retries: int = 3,
    mode: str = "local",
    run_id: str | None = None,
    synth_mode: str = "kl",
) -> dict:
    run_id = run_id or new_run_id(taxon_slug)
    history: list[dict] = []
    feedback: list | str | None = None

    for attempt in range(1, max_retries + 1):
        hyp = write_hypothesis_bundle(taxon_slug, run_id, retry_feedback=feedback)
        history.append({"attempt": attempt, "state": "HYPOTHESIS", "prompt": hyp["prompt_path"]})

        if mode == "agents":
            result = {
                "ok": False,
                "mode": mode,
                "state": "AWAITING_ANTIGRAVITY",
                "run_id": run_id,
                "taxon": taxon_slug,
                "attempt": attempt,
                "antigravity_prompt": hyp["prompt_path"],
                "history": history,
                "instruction": (
                    "Delegate BUILD_AND_SIM to Antigravity with the prompt file, then call "
                    f"`python -m pipeline.orchestrator --taxon {taxon_slug} --run-id {run_id} --resume-review`."
                ),
            }
            _write_summary(run_id, result)
            return result

        metrics = run_build_and_sim(taxon_slug, run_id, synth_mode=synth_mode)
        history.append({"attempt": attempt, "state": "BUILD_AND_SIM", "metrics_state": metrics.get("state")})

        review = review_run(run_id)
        term_state = normalize_legacy_state(review["terminal_state"])
        history.append(
            {
                "attempt": attempt,
                "state": "PEER_REVIEW",
                "passed": review["passed"],
                "terminal_state": term_state,
                "issues": review["issues"],
            }
        )

        if review["passed"]:
            result = {
                "ok": True,
                "mode": mode,
                "state": term_state,
                "contract_status": review["contract_status"],
                "run_id": run_id,
                "taxon": taxon_slug,
                "attempt": attempt,
                "output_dir": str(OUTPUT / run_id),
                "history": history,
                "review": review,
            }
            _write_summary(run_id, result)
            return result

        feedback = review.get("retry_feedback")
        history.append({"attempt": attempt, "state": "RETRY_LOOP", "feedback": feedback})

    result = {
        "ok": False,
        "mode": mode,
        "state": "REVIEW_FAILED",
        "run_id": run_id,
        "taxon": taxon_slug,
        "attempt": max_retries,
        "output_dir": str(OUTPUT / run_id),
        "history": history,
        "last_feedback": feedback,
    }
    _write_summary(run_id, result)
    return result


def resume_review(taxon_slug: str, run_id: str, *, max_retries: int = 3, attempt: int = 1) -> dict:
    """After Antigravity finished BUILD_AND_SIM, Grokbot/local peer-review continues."""
    review = review_run(run_id)
    term_state = normalize_legacy_state(review["terminal_state"])
    if review["passed"]:
        result = {
            "ok": True,
            "state": term_state,
            "contract_status": review["contract_status"],
            "run_id": run_id,
            "taxon": taxon_slug,
            "review": review,
            "output_dir": str(OUTPUT / run_id),
        }
        _write_summary(run_id, result)
        return result

    if attempt >= max_retries:
        result = {
            "ok": False,
            "state": "REVIEW_FAILED",
            "run_id": run_id,
            "taxon": taxon_slug,
            "review": review,
            "last_feedback": review.get("retry_feedback"),
        }
        _write_summary(run_id, result)
        return result

    hyp = write_hypothesis_bundle(taxon_slug, run_id, retry_feedback=review.get("retry_feedback"))
    result = {
        "ok": False,
        "state": "AWAITING_ANTIGRAVITY_RETRY",
        "run_id": run_id,
        "taxon": taxon_slug,
        "attempt": attempt + 1,
        "antigravity_prompt": hyp["prompt_path"],
        "review": review,
    }
    _write_summary(run_id, result)
    return result


def _write_summary(run_id: str, result: dict) -> None:
    path = OUTPUT / run_id / "pipeline_summary.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Paleo HYPOTHESIS->BUILD_AND_SIM->PEER_REVIEW orchestrator")
    parser.add_argument("--taxon", default="parasaurolophus_ram14000")
    parser.add_argument("--list-taxa", action="store_true")
    parser.add_argument("--mode", choices=["local", "agents"], default="local")
    parser.add_argument("--synth", choices=["kl", "legacy_iir"], default="kl")
    parser.add_argument("--max-retries", type=int, default=3)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--resume-review", action="store_true")
    parser.add_argument("--attempt", type=int, default=1)
    args = parser.parse_args(argv)

    if args.list_taxa:
        print(json.dumps(list_taxa(), indent=2))
        return

    if args.resume_review:
        if not args.run_id:
            raise SystemExit("--run-id is required with --resume-review")
        result = resume_review(args.taxon, args.run_id, max_retries=args.max_retries, attempt=args.attempt)
    else:
        result = run_pipeline(
            args.taxon,
            max_retries=args.max_retries,
            mode=args.mode,
            run_id=args.run_id,
            synth_mode=args.synth,
        )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
