find pvol | awk '{filein=$0;gsub(/pvol/,"vp")}{system("vol2bird "filein" "$0)}'
