cat -n plus_minus_report.map | sort -uk2 | sort -nk1 | cut -f2- > plus_minus_report_unique.map
mv plus_minus_report_unique.map plus_minus_report.map
