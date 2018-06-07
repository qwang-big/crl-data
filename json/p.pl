sub x {
return '{"type":"bigwig","name":"'.$_[0].'","url":"'.$_[1].'","mode":1,"qtc":{"anglescale":1,"pr":0,"pg":0,"pb":230,"nr":255,"ng":0,"nb":0,"pth":"#000099","nth":"#800000","thtype":0,"thmin":0,"thmax":10,"thpercentile":90,"height":20,"summeth":1},"metadata":{}},'."\n"
}
while(<>){
$s=$1 if /"big_data_url": "(\S+)"/i;
$m=$1 if /"assay": "(\S+)",/i;
$h{$m}.=x("H1-$m",$s) if /"cell_type": "H1",/i;
$h{$m}.=x("NPC-$m",$s) if /"cell_type": "H1_Derived_Neuronal_Progenitor_Cultured_Cells",/i;
#$h{$m}.=x("MSC-$m",$s) if /"CELL_TYPE": "H1_Derived_Mesenchymal_Stem_Cells",/i;
#$h{$m}.=x("TSC-$m",$s) if /"CELL_TYPE": "H1_BMP4_Derived_Trophoblast_Cultured_Cells",/i;
}
print '[';
foreach(sort keys %h){
print $h{$_}
}
print '{"type":"native_track","list":[{"name":"refGene","mode":3,"qtc":{"anglescale":1,"pr":255,"pg":0,"pb":0,"nr":0,"ng":0,"nb":230,"pth":"#800000","nth":"#000099","thtype":0,"thmin":0,"thmax":10,"thpercentile":95,"height":50,"summeth":1,"textcolor":"#000000","fontsize":"8pt","fontfamily":"sans-serif","fontbold":false,"bedcolor":"#002EB8"},"metadata":{},"details":{"source":"UCSC Genome Browser","download date":"Jan 1, 2012"}}]},{"type":"metadata","vocabulary_set":{}}]'
