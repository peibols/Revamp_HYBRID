#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest


MODULE_PATH = Path(__file__).with_name("convert_oo_paired_to_root.py")
SPEC = importlib.util.spec_from_file_location("convert_oo_paired_to_root", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
converter = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = converter
SPEC.loader.exec_module(converter)


def event_text(px: float, weight: float = 2.5) -> bytes:
    return (
        "# event 0\n"
        f"weight {weight} cross 17.5 X 0.2 Y -0.3\n"
        "5 0 0 0 21 -2\n"
        f"{px} 0 0 0.13957 211 0\n"
        "1.0 0.2 0 0.13957 211 1\n"
        "0.1 0.02 0 0.13957 -211 2\n"
        "end\n"
    ).encode()


def add_member(tar: tarfile.TarFile, name: str, payload: bytes) -> None:
    member = tarfile.TarInfo(name)
    member.size = len(payload)
    tar.addfile(member, io.BytesIO(payload))


def make_archive(
    path: Path,
    with_weight: float = 2.5,
    *,
    task_id: int = 7,
    seed: int = 12345,
    hydro_index: int = 3,
) -> None:
    base = "runs/campaign/aa/hydro_03_C0-5"
    payload_sha256 = "a" * 64
    task = f"task_{task_id:05d}"
    pair_summary = (
        "kind\ttask_id\tseed\tevents\tcentrality\thydro_index\thydro_event_id"
        "\thydro_ncoll\thydro_payload_sha256\tvariant\tuse_prehydro"
        "\tprehydro_file\treturncode\ttimeout\tseconds\tdir\n"
        f"aa\t{task_id}\t{seed}\t1\tC0-5\t{hydro_index}\t777\t42\t{payload_sha256}"
        "\tno_prehydro\t0\t\t0\t0\t1.0\t/no\n"
        f"aa\t{task_id}\t{seed}\t1\tC0-5\t{hydro_index}\t777\t42\t{payload_sha256}"
        "\twith_prehydro\t1\tpre.tsv\t0\t0\t1.1\t/with\n"
    ).encode()
    summary_header = "variant\tuse_prehydro\tprehydro_file\treturncode\ttimeout\tseconds\tdir\n"
    with tarfile.open(path, "w:gz") as tar:
        add_member(tar, f"{base}/{task}_pair_summary.tsv", pair_summary)
        add_member(
            tar,
            f"{base}/{task}/summary.tsv",
            (summary_header + "no_prehydro\t0\t\t0\t0\t1.0\t/no\n").encode(),
        )
        add_member(tar, f"{base}/{task}/HYBRID_Hadrons.out", event_text(4.0))
        add_member(
            tar,
            f"{base}/{task}_prehydro/summary.tsv",
            (summary_header + "with_prehydro\t1\tpre.tsv\t0\t0\t1.1\t/with\n").encode(),
        )
        add_member(
            tar,
            f"{base}/{task}_prehydro/HYBRID_Hadrons.out",
            event_text(3.5, with_weight),
        )


class ConverterTest(unittest.TestCase):
    def test_parse_strict_pair(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "chunk_7.tar.gz"
            make_archive(archive)
            pair = converter.parse_paired_archive(archive, expected_chunk_id=7)
        self.assertEqual(pair.seed, 12345)
        self.assertEqual(pair.hydro_index, 3)
        self.assertEqual(pair.hydro_event_id, 777)
        self.assertEqual(pair.hydro_ncoll, 42)
        self.assertEqual(pair.hydro_payload_sha256, "a" * 64)
        self.assertEqual(pair.no_prehydro.event_number, 0)
        self.assertEqual(len(pair.no_prehydro.particles), 4)
        self.assertEqual([particle.raw_label for particle in pair.no_prehydro.particles], [-2, 0, 1, 2])

    def test_reject_variant_metadata_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "chunk_7.tar.gz"
            make_archive(archive, with_weight=3.0)
            with self.assertRaisesRegex(converter.ArchiveValidationError, "event weight"):
                converter.parse_paired_archive(archive, expected_chunk_id=7)

    def test_pair_id_distinguishes_overlapping_sources(self) -> None:
        self.assertEqual(converter.stable_pair_id(0, 7), 7)
        self.assertEqual(converter.stable_pair_id(1, 7), (1 << 48) | 7)
        self.assertNotEqual(converter.stable_pair_id(0, 7), converter.stable_pair_id(1, 7))

    def test_validate_pair_against_v2_task_assignment(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = Path(temporary) / "chunk_7.tar.gz"
            make_archive(archive)
            pair = converter.parse_paired_archive(archive, expected_chunk_id=7)
        assignment = {
            "hard_seed": "12345",
            "hydro_slot": "3",
            "hydro_event_id": "777",
            "hydro_ncoll": "42",
            "hydro_payload_sha256": "a" * 64,
        }
        converter.validate_pair_task_assignment(
            pair, task_id=7, assignment=assignment
        )
        assignment["hydro_event_id"] = "778"
        with self.assertRaisesRegex(
            converter.ArchiveValidationError, "hydro_event_id"
        ):
            converter.validate_pair_task_assignment(
                pair, task_id=7, assignment=assignment
            )

    @unittest.skipUnless(
        shutil.which("root-config") and shutil.which("fastjet-config"),
        "ROOT and FastJet are required",
    )
    def test_end_to_end_root_writer(self) -> None:
        import awkward as ak
        import uproot

        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "source"
            (source / "outputs" / "aa").mkdir(parents=True)
            (source / "status" / "aa").mkdir(parents=True)
            make_archive(source / "outputs" / "aa" / "chunk_7.tar.gz")
            (source / "status" / "aa" / "chunk_7.txt").write_text(
                "kind=aa\ntask_id=7\nstatus=success\nexit_code=0\n"
            )
            output = work / "paired.root"
            manifest = work / "aa_task_manifest.tsv"
            manifest.write_text(
                "task_id\thard_seed\thydro_slot\thydro_event_id\thydro_ncoll"
                "\thydro_payload_sha256\n"
                f"7\t12345\t3\t777\t42\t{'a' * 64}\n"
            )
            subprocess.run(
                [
                    sys.executable,
                    str(MODULE_PATH),
                    "--source",
                    f"synthetic={source}",
                    "--output",
                    str(output),
                    "--build-dir",
                    str(work / "build"),
                    "--aa-task-manifest",
                    str(manifest),
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            with uproot.open(output) as root_file:
                hadrons = root_file["noPrehydro/Hadrons"].arrays(library="ak")
                self.assertEqual(ak.to_list(hadrons.hadronStatus[0]), [0, 1, -1])
                self.assertEqual(ak.to_list(hadrons.hadronRawLabel[0]), [0, 1, 2])
                self.assertEqual(int(hadrons.nHardMarkers[0]), 1)
                jets = root_file["noPrehydro/Jets"].arrays(
                    [
                        "jet1Pt",
                        "jet1RawPt",
                        "jet1NegativeWakePt",
                        "jet2Pt",
                        "jet2RawPt",
                        "jet2NegativeWakePt",
                        "jet2NormalPtD",
                        "jet2NormalEffectiveMultiplicity",
                        "jet2LeadingNormalFraction",
                        "jet2HardPartonId",
                        "jet2HardPartonPt",
                        "jet2HardPartonDR",
                        "jet2PairMatchIndex",
                        "jet2PairMatchDR",
                        "jet2PairMatchOtherPt",
                        "jet2PairMatchOtherHardPartonId",
                        "jet4Pt",
                        "jet4RawPt",
                        "jet4NegativeWakePt",
                    ],
                    library="ak",
                )
                self.assertGreater(float(jets.jet1RawPt[0][0]), 0.0)
                self.assertGreaterEqual(float(jets.jet1NegativeWakePt[0][0]), 0.0)
                self.assertLessEqual(float(jets.jet1Pt[0][0]), float(jets.jet1RawPt[0][0]))
                self.assertGreater(float(jets.jet2NegativeWakePt[0][0]), 0.0)
                self.assertLess(float(jets.jet2Pt[0][0]), float(jets.jet2RawPt[0][0]))
                self.assertAlmostEqual(float(jets.jet2NormalPtD[0][0]), 1.0)
                self.assertAlmostEqual(
                    float(jets.jet2NormalEffectiveMultiplicity[0][0]), 1.0
                )
                self.assertAlmostEqual(float(jets.jet2LeadingNormalFraction[0][0]), 1.0)
                self.assertEqual(int(jets.jet2HardPartonId[0][0]), 21)
                self.assertAlmostEqual(float(jets.jet2HardPartonPt[0][0]), 5.0)
                self.assertLess(float(jets.jet2HardPartonDR[0][0]), 0.2)
                self.assertEqual(int(jets.jet2PairMatchIndex[0][0]), 0)
                self.assertLess(float(jets.jet2PairMatchDR[0][0]), 0.1)
                self.assertGreater(float(jets.jet2PairMatchOtherPt[0][0]), 0.0)
                self.assertEqual(int(jets.jet2PairMatchOtherHardPartonId[0][0]), 21)
                self.assertGreater(float(jets.jet4NegativeWakePt[0][0]), 0.0)
                self.assertLess(float(jets.jet4Pt[0][0]), float(jets.jet4RawPt[0][0]))

    @unittest.skipUnless(
        shutil.which("root-config") and shutil.which("fastjet-config"),
        "ROOT and FastJet are required",
    )
    def test_combined_manifest_validates_two_disjoint_sources(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            sources = []
            for name, task_id, seed, hydro_index in (
                ("first", 7, 12345, 3),
                ("continuation", 8, 12346, 4),
            ):
                source = work / name
                (source / "outputs/aa").mkdir(parents=True)
                (source / "status/aa").mkdir(parents=True)
                make_archive(
                    source / f"outputs/aa/chunk_{task_id}.tar.gz",
                    task_id=task_id,
                    seed=seed,
                    hydro_index=hydro_index,
                )
                (source / f"status/aa/chunk_{task_id}.txt").write_text(
                    f"kind=aa\ntask_id={task_id}\nstatus=success\nexit_code=0\n"
                )
                sources.append(source)

            manifest = work / "aa_task_manifest.tsv"
            manifest.write_text(
                "task_id\thard_seed\thydro_slot\thydro_event_id\thydro_ncoll"
                "\thydro_payload_sha256\n"
                f"7\t12345\t3\t777\t42\t{'a' * 64}\n"
                f"8\t12346\t4\t777\t42\t{'a' * 64}\n"
            )
            output = work / "paired.root"
            subprocess.run(
                [
                    sys.executable,
                    str(MODULE_PATH),
                    "--source",
                    f"first={sources[0]}",
                    "--source",
                    f"continuation={sources[1]}",
                    "--output",
                    str(output),
                    "--build-dir",
                    str(work / "build"),
                    "--aa-task-manifest",
                    str(manifest),
                ],
                check=True,
                text=True,
                capture_output=True,
            )
            summary = json.loads(output.with_suffix(".summary.json").read_text())
            self.assertEqual(summary["acceptedPairs"], 2)
            self.assertEqual(summary["aaTaskManifest"]["acceptedRows"], 2)
            self.assertEqual(summary["aaTaskManifest"]["closure"], "PASS")

    @unittest.skipUnless(
        shutil.which("root-config") and shutil.which("fastjet-config"),
        "ROOT and FastJet are required",
    )
    def test_manifest_closure_rejects_missing_task(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            source = work / "source"
            (source / "outputs/aa").mkdir(parents=True)
            (source / "status/aa").mkdir(parents=True)
            make_archive(source / "outputs/aa/chunk_7.tar.gz")
            (source / "status/aa/chunk_7.txt").write_text(
                "kind=aa\ntask_id=7\nstatus=success\nexit_code=0\n"
            )
            manifest = work / "aa_task_manifest.tsv"
            manifest.write_text(
                "task_id\thard_seed\thydro_slot\thydro_event_id\thydro_ncoll"
                "\thydro_payload_sha256\n"
                f"7\t12345\t3\t777\t42\t{'a' * 64}\n"
                f"8\t12346\t4\t777\t42\t{'a' * 64}\n"
            )
            output = work / "paired.root"
            result = subprocess.run(
                [
                    sys.executable,
                    str(MODULE_PATH),
                    "--source",
                    f"first={source}",
                    "--output",
                    str(output),
                    "--build-dir",
                    str(work / "build"),
                    "--aa-task-manifest",
                    str(manifest),
                ],
                check=False,
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("AA task-manifest closure failed", result.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
