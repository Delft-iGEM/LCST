"""Command-line interface for ELP LCST prediction."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .model import predict
from .sequence import guests_from_repeat_spec


def _spec_to_sequence(guest_spec: str, n_repeats: int) -> str:
    guests = guests_from_repeat_spec(guest_spec) * n_repeats
    return "".join(f"VPG{g}G" for g in guests)


def _cmd_predict(args) -> int:
    if args.spec:
        if args.repeats is None:
            print("error: --repeats is required with --spec", file=sys.stderr)
            return 2
        sequence = _spec_to_sequence(args.spec, args.repeats)
    elif args.sequence:
        sequence = args.sequence
    else:
        sequence = sys.stdin.read().strip()

    pred = predict(sequence, concentration_uM=args.concentration, pH=args.pH)

    if args.json:
        print(json.dumps(pred.as_dict(), indent=2))
        return 0

    print(f"Predicted LCST (transition temperature):  {pred.tt_celsius:.1f} degC")
    print(f"  chain length            : {pred.n_pentads} pentapeptides")
    print(f"  molecular weight        : {pred.molecular_weight_da/1000:.1f} kDa")
    print(f"  mean guest Urry T_tc    : {pred.mean_guest_ttc:.1f} degC")
    print(f"  effective Ala fraction  : {pred.effective_fraction_ala:.2f}")
    print(f"  fraction charged guests : {pred.fraction_charged:.2f}")
    print(f"  fraction aromatic guests: {pred.fraction_aromatic:.2f}")
    print(f"  concentration           : {pred.concentration_uM:g} uM")
    print(f"  pH                      : {pred.pH:g}")
    print(f"  recognised as ELP       : {pred.is_elp} "
          f"({pred.motif_match_fraction:.0%} pentads match VPGXG)")
    for w in pred.warnings:
        print(f"  ! {w}")
    return 0


def _cmd_validate(args) -> int:
    from .validate import report
    print(report())
    return 0


def _cmd_calibrate(args) -> int:
    import pandas as pd
    from . import calibrate

    df_raw = pd.read_csv(args.data)
    if "sequence" in df_raw.columns:
        train_df = calibrate.dataset_from_sequences(df_raw)
    else:
        train_df = df_raw
    if args.with_unified:
        train_df = pd.concat(
            [calibrate.generate_unified_dataset(), train_df], ignore_index=True
        )
    cm = calibrate.train(train_df)
    cm.save(args.out)
    print(f"Trained on {cm.n_train} rows.")
    print(f"  cross-validated MAE : {cm.cv_mae:.2f} degC")
    print(f"  cross-validated RMSE: {cm.cv_rmse:.2f} degC")
    print(f"Saved calibrated model to {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="elp-lcst",
        description="Predict the LCST (inverse transition temperature) of "
                    "elastin-like polypeptide (ELP) sequences.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    pr = sub.add_parser("predict", help="predict Tt for a sequence")
    g = pr.add_mutually_exclusive_group()
    g.add_argument("sequence", nargs="?", help="one-letter amino-acid sequence")
    g.add_argument("--spec", help="compact guest spec, e.g. V5A2G3")
    pr.add_argument("--repeats", type=int, help="repeats of --spec block")
    pr.add_argument("-c", "--concentration", type=float, default=25.0,
                    help="ELP concentration in uM (default 25)")
    pr.add_argument("--pH", type=float, default=7.4, help="solution pH (default 7.4)")
    pr.add_argument("--json", action="store_true", help="emit JSON")
    pr.set_defaults(func=_cmd_predict)

    pv = sub.add_parser("validate", help="report model accuracy vs literature")
    pv.set_defaults(func=_cmd_validate)

    pc = sub.add_parser("calibrate", help="train ML residual correction from data")
    pc.add_argument("data", help="CSV with measured Tt (see README for columns)")
    pc.add_argument("-o", "--out", default="calibrated_model.pkl",
                    help="output model path")
    pc.add_argument("--with-unified", action="store_true",
                    help="augment with the unified-model V/A grid")
    pc.set_defaults(func=_cmd_calibrate)

    args = p.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
