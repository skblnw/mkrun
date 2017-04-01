package require Orient
namespace import Orient::orient

set sel [atomselect top "all"]

set A [transaxis y 180]
$sel move $A