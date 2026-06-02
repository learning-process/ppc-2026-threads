from main import _has_single_nonempty_markdown, _has_thread_report


def test_has_single_nonempty_markdown_requires_exactly_one_file(temp_dir):
    report_dir = temp_dir / "task"
    report_dir.mkdir()

    assert not _has_single_nonempty_markdown(report_dir)

    report_path = report_dir / "report.md"
    report_path.write_text("\n")
    assert not _has_single_nonempty_markdown(report_dir)

    report_path.write_text("report\n")
    assert _has_single_nonempty_markdown(report_dir)

    (report_dir / "extra.md").write_text("extra\n")
    assert not _has_single_nonempty_markdown(report_dir)


def test_has_thread_report_requires_root_and_implementation_reports(temp_dir):
    task_dir = temp_dir / "task"
    impl_dir = task_dir / "omp"
    impl_dir.mkdir(parents=True)

    (task_dir / "report.md").write_text("root report\n")
    assert not _has_thread_report(task_dir, "omp")

    (impl_dir / "report.md").write_text("omp report\n")
    assert _has_thread_report(task_dir, "omp")

    (impl_dir / "extra.md").write_text("extra\n")
    assert not _has_thread_report(task_dir, "omp")
