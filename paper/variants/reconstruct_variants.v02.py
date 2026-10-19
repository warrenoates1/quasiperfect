"""Reconstruct every version of every source file from the session transcript.

Why this exists: several intermediate versions of the solver were written,
tested, found wrong and corrected before ever being committed, so they are in
no git commit. Claude Code's own file-history backups do not hold them either
(the backup filenames are null for files created during the session). What does
survive is the transcript, which records the full text of every Write and the
exact old/new strings of every Edit.

Replaying those in order reconstructs each version exactly. The replay is
self-checking: if the final reconstructed state of a file matches the file now
on disk byte for byte, every earlier state in that chain is also faithful,
because each step is an exact textual operation. Any file whose final state does
not match is reported as UNVERIFIED and must not be trusted.

Usage: python reconstruct_variants.py SESSION.jsonl [...]   (chronological order)
Writes variants/<name>.vNN<suffix> and variants/MANIFEST.md.
"""
import hashlib
import json
import os
import sys

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
OUT = os.path.join(os.path.dirname(__file__), "variants")

# Short labels tying a version to the issue it represents in the paper. These are
# filled in only after reading the reconstructed diffs: guessing which version is
# which would put an unchecked claim into the record. Keyed by (file, version).
LABELS = {}


def rel(p):
    p = p.replace("/", os.sep)
    return os.path.relpath(os.path.abspath(p), REPO) if os.path.isabs(p) else p


def replay(paths):
    state, versions = {}, {}
    for path in paths:
        for line in open(path, encoding="utf-8"):
            d = json.loads(line)
            if d.get("type") != "assistant":
                continue
            ts = d.get("timestamp", "")
            for b in d["message"].get("content", []):
                if b.get("type") != "tool_use":
                    continue
                name, inp = b["name"], b.get("input", {})
                if name not in ("Write", "Edit") or "file_path" not in inp:
                    continue
                f = rel(inp["file_path"])
                if not f.endswith(".py"):
                    continue
                if name == "Write":
                    state[f] = inp.get("content", "")
                    op = "Write (full file)"
                else:
                    old, new = inp.get("old_string", ""), inp.get("new_string", "")
                    if f not in state or old not in state[f]:
                        versions.setdefault(f, []).append(
                            (ts, "Edit FAILED TO REPLAY", None))
                        continue
                    n = -1 if inp.get("replace_all") else 1
                    state[f] = state[f].replace(old, new, n) if n == 1 else \
                        state[f].replace(old, new)
                    op = "Edit"
                versions.setdefault(f, []).append((ts, op, state[f]))
    return state, versions


def main():
    state, versions = replay(sys.argv[1:])
    os.makedirs(OUT, exist_ok=True)
    lines = ["# Reconstructed source variants", "",
             "Each version below is the exact state of the file after one recorded",
             "Write or Edit, replayed from the session transcript. `verified` means the",
             "final replayed state matches the file currently on disk byte for byte,",
             "which makes the whole chain for that file faithful.", ""]
    for f in sorted(versions):
        base = os.path.basename(f)
        stem, ext = os.path.splitext(base)
        disk = os.path.join(REPO, f)
        final = state.get(f)
        ok = os.path.exists(disk) and open(disk, encoding="utf-8").read() == final
        lines += [f"## {f} — {'verified' if ok else '**UNVERIFIED**'}", "",
                  "| v | time (UTC) | operation | sha1 | label |",
                  "|---|---|---|---|---|"]
        for i, (ts, op, text) in enumerate(versions[f], 1):
            if text is None:
                lines.append(f"| {i} | {ts[:19]} | {op} | – | – |")
                continue
            name = f"{stem}.v{i:02d}{ext}"
            open(os.path.join(OUT, name), "w", encoding="utf-8", newline="\n").write(text)
            sha = hashlib.sha1(text.encode()).hexdigest()[:12]
            lines.append(f"| {i} | {ts[:19]} | {op} | `{sha}` | "
                         f"{LABELS.get((base, i), '')} |")
        lines.append("")
        print(f"{f}: {len(versions[f])} versions, "
              f"{'verified' if ok else 'UNVERIFIED — do not trust'}")
    open(os.path.join(OUT, "MANIFEST.md"), "w", encoding="utf-8").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
