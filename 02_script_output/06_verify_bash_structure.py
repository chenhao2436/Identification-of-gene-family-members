#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
轻量 bash 结构校验器（替代无法在沙箱内运行的 `bash -n`）。

背景：本机 Windows 沙箱禁止 Git Bash（无法创建 signal pipe），
      `bash -n` 无法执行，故用静态检查兜底。

检查项
    1. () {} 配对（跳过单引号、双引号、注释、heredoc 内部）
    2. do/done、if/fi、case/esac 关键字配对
    3. heredoc 起始与结束标记配对
    4. 行尾统一为 LF（CRLF 会让 Linux 上 bash 报错）
    5. 危险模式：`$(... )` 内未加保护的 grep -c（set -e 下会提前中止）
"""
import os
import re
import sys

KEYWORDS = [("if", "fi"), ("do", "done"), ("case", "esac")]

# case 分支模式标签：`--query)`、`*)`、`-h|--help)` 等，其中的 ')' 不是语法错误。
# 注意排除 `;;`（分支结束）、`$`（命令替换）、花括号（函数体），这些不是标签。
CASE_LABEL_RE = re.compile(r"^\s*(?:[^\s()|]+\|)*[^\s()|]+\)(?!\s*;)")


def is_case_label(code):
    """判断剥离引号后的代码行是否为 case 分支标签"""
    s = code.strip()
    if not s or s.startswith("$") or s.startswith("{") or s.startswith("}"):
        return False
    return bool(CASE_LABEL_RE.match(s))


def strip_quotes_stateful(line, in_single):
    """
    剥离单引号/双引号内容；单引号状态跨行保持（bash 允许跨行单引号串）。
    返回 (剥离后的代码, 行末是否仍在单引号内)。
    """
    out = []
    i = 0
    in_d = False
    n = len(line)
    while i < n:
        ch = line[i]
        if in_single:
            if ch == "'":
                in_single = False
            i += 1
            continue
        if ch == "'" and not in_d:
            in_single = True
            i += 1
            continue
        if ch == '"' and not in_single:
            in_d = not in_d
            i += 1
            continue
        if ch == "\\" and i + 1 < n:
            i += 2
            continue
        # bash 里 `#` 只有在“词首”（前一个字符是空白或行首）才是注释起始。
        # 否则 `$#`（参数个数）、`${VAR#pat}` 会被误当注释，截断代码。
        if ch == "#" and not in_d and (i == 0 or line[i - 1] in " \t;"):
            break
        out.append(ch)
        i += 1
    return "".join(out), in_single


def check_file(path):
    with open(path, "rb") as fh:
        raw = fh.read()

    problems = []
    notes = []

    # ---- 行尾检查 ----
    crlf = raw.count(b"\r\n")
    if crlf:
        problems.append("含 %d 个 CRLF 行尾（Linux 上 bash 会报错），需转 LF" % crlf)
    if b"\r" in raw.replace(b"\r\n", b""):
        problems.append("含孤立 CR 字符")

    text = raw.decode("utf-8", errors="replace")
    lines = text.split("\n")

    # ---- 逐行扫描，剥离注释/引号/heredoc ----
    depth_paren = 0
    depth_brace = 0
    kw_stack = []
    heredoc = None
    heredoc_count = 0
    in_single = False          # 跨行单引号状态

    for ln, line in enumerate(lines, 1):
        stripped = line.strip()

        # heredoc 内部原样跳过
        if heredoc is not None:
            if stripped == heredoc:
                heredoc = None
                heredoc_count += 1
            continue

        # 处于跨行单引号串内部（如 usage() 的文档块），整行跳过
        if in_single:
            _, in_single = strip_quotes_stateful(line, True)
            continue

        # 检测 heredoc 起始（<<EOF / <<'EOF' / <<"EOF"）
        m = re.search(r"<<-?\s*['\"]?([A-Za-z_][A-Za-z0-9_]*)['\"]?", line)
        if m and "<<" in line and not line.lstrip().startswith("#"):
            heredoc = m.group(1)
            continue

        code, in_single = strip_quotes_stateful(line, False)

        # case 分支标签行：`--query)` / `*)` 里的 ')' 是标签结束符，不是语法闭合；
        # 标签后的消息文本里也可能含括号（如 `die "...(...)..."`），一并跳过，
        # 否则会造成假的不配对。
        case_label = is_case_label(code)

        if not case_label:
            depth_paren += code.count("(") - code.count(")")
            if depth_paren < 0:
                problems.append("第 %d 行: 右括号 ')' 多于 '('" % ln)
                depth_paren = 0

        depth_brace += code.count("{") - code.count("}")
        if depth_brace < 0:
            problems.append("第 %d 行: 右花括号 '}' 多于 '{'" % ln)
            depth_brace = 0

        # 关键字配对（只看剥离后的代码 token）
        toks = re.findall(r"\b(if|fi|do|done|case|esac|for|while|until)\b", code)
        # `while ...` / `for ...` 与 `do` 分处两行时，块开始是那行 `do`，
        # 此时不要把循环关键字当作块开始，否则 done 会配不上。
        if ("for" in toks or "while" in toks or "until" in toks) and "do" not in toks:
            toks = [t for t in toks if t not in ("for", "while", "until")]

        for t in toks:
            if t in ("if", "do", "case"):
                kw_stack.append((t, ln))
            elif t == "fi":
                if kw_stack and kw_stack[-1][0] == "if":
                    kw_stack.pop()
                else:
                    problems.append("第 %d 行: 多余的 'fi'（无匹配 if）" % ln)
            elif t == "done":
                if kw_stack and kw_stack[-1][0] == "do":
                    kw_stack.pop()
                else:
                    problems.append("第 %d 行: 多余的 'done'（无匹配 do）" % ln)
            elif t == "esac":
                if kw_stack and kw_stack[-1][0] == "case":
                    kw_stack.pop()
                else:
                    problems.append("第 %d 行: 多余的 'esac'（无匹配 case）" % ln)

    if heredoc is not None:
        problems.append("heredoc 未闭合，缺少结束标记 '%s'" % heredoc)
    if depth_paren != 0:
        problems.append("圆括号不配对，最终深度 %d" % depth_paren)
    if depth_brace != 0:
        problems.append("花括号不配对，最终深度 %d" % depth_brace)
    for kw, ln in kw_stack:
        problems.append("第 %d 行的 '%s' 没有对应的结束关键字" % (ln, kw))

    # ---- set -e 危险模式检查 ----
    for ln, line in enumerate(lines, 1):
        if line.lstrip().startswith("#"):
            continue
        # $( ... grep -c ... ) 未带 || true
        for m in re.finditer(r"\$\(([^)]*grep\s+-c[^)]*)\)", line):
            inner = m.group(1)
            if "|| true" not in inner and "||true" not in inner:
                problems.append(
                    "第 %d 行: $(%s) 未加 '|| true'，计数为 0 时 grep 返回 1，"
                    "在 set -e 下会提前中止脚本" % (ln, inner.strip()))

    notes.append("总行数: %d" % len(lines))
    notes.append("heredoc 块: %d 个" % heredoc_count)
    return problems, notes


def check_defaults(script_path):
    """
    复刻流水线的「零参数默认路径」探测逻辑（BASH_SOURCE -> SCRIPT_DIR -> PROJECT_DIR），
    并核对磁盘上这些路径是否真的存在、序列条数是否正确。

    沙箱内无法运行 bash，因此用本函数验证 --query/--db-fasta 的默认值是否可用，
    避免"一键运行"实际上找不到文件。
    """
    problems = []
    script_path = os.path.abspath(script_path)
    script_dir = os.path.dirname(script_path)
    project_dir = os.path.dirname(script_dir)

    print("=" * 66)
    print("默认路径检查（模拟零参数运行）")
    print("=" * 66)
    print("  SCRIPT_DIR  = %s" % script_dir)
    print("  PROJECT_DIR = %s" % project_dir)

    # 与脚本中 probe 顺序一致
    query_cands = [
        os.path.join(script_dir, "01_At_TRM_query_34.fasta"),
        os.path.join(project_dir, "01_At_TRM_query_34.fasta"),
        os.path.join(project_dir, "01_At_TRM_query", "01_At_TRM_query_34.fasta"),
    ]
    input_dir = os.path.join(project_dir, "01_At_TRM_query")
    db_cands = [
        os.path.join(input_dir, "Zunla_Canz.pep.fa"),
        os.path.join(input_dir, "Canz.pep.fa"),
        os.path.join(project_dir, "Zunla_Canz.pep.fa"),
    ]
    map_cands = [
        os.path.join(script_dir, "02_At_TRM_annotation.tsv"),
        os.path.join(project_dir, "02_At_TRM_annotation.tsv"),
    ]
    outdir = os.path.join(project_dir, "05_blast_results")

    def first(cands, label, expect_headers=None):
        for c in cands:
            if os.path.isfile(c):
                n = None
                if expect_headers is not None:
                    with open(c, "r", encoding="utf-8", errors="replace") as fh:
                        n = sum(1 for l in fh if l.startswith(">"))
                    ok = (n == expect_headers)
                    mark = "OK" if ok else "条数不符!"
                    print("  [%s] %-10s %s  (%d 条, 期望 %d)"
                          % (mark, label, c, n, expect_headers))
                    if not ok:
                        problems.append("%s 序列数为 %d，期望 %d（%s）"
                                        % (label, n, expect_headers, c))
                else:
                    print("  [OK] %-10s %s" % (label, c))
                return c
        problems.append("%s 的默认路径全部不存在，尝试过：%s"
                        % (label, "; ".join(cands)))
        print("  [FAIL] %-10s 未找到（尝试 %d 个候选路径）" % (label, len(cands)))
        return None

    first(query_cands, "query", expect_headers=34)
    first(db_cands, "db", expect_headers=52385)
    first(map_cands, "map")
    print("  [--] outdir     %s  (脚本会自动创建)" % outdir)

    print("-" * 66)
    if problems:
        print("默认路径检查失败 %d 项:" % len(problems))
        for p in problems:
            print("  - %s" % p)
        return problems
    print("默认路径检查通过 ✅  可直接零参数一键运行")
    return []


def main():
    argv = sys.argv[1:]
    if not argv:
        print("用法: python verify_bash_structure.py <script.sh> [...]")
        print("      python verify_bash_structure.py --check-defaults <pipeline.sh>")
        return 2

    if argv[0] == "--check-defaults":
        if len(argv) < 2:
            print("用法: python verify_bash_structure.py --check-defaults <pipeline.sh>")
            return 2
        if not os.path.exists(argv[1]):
            print("[错误] 找不到文件: %s" % argv[1])
            return 1
        return 1 if check_defaults(argv[1]) else 0

    rc = 0
    for path in argv:
        if not os.path.exists(path):
            print("[错误] 找不到文件: %s" % path)
            rc = 1
            continue
        problems, notes = check_file(path)
        print("=" * 66)
        print("检查: %s" % path)
        for n in notes:
            print("  %s" % n)
        if problems:
            rc = 1
            print("  发现问题 %d 个:" % len(problems))
            for p in problems:
                print("    [FAIL] %s" % p)
        else:
            print("  结构检查通过 ✅")
    return rc


if __name__ == "__main__":
    sys.exit(main())
