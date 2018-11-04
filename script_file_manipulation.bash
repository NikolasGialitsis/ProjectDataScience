printf '\n%.0s' {0..999} >> filename

awk '{ print $0NR }' filename >> filename2

paste -d' ' filename2 tiny.dat >> final_dataset