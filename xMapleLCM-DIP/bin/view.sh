# This is for displaying the conflict graphs that we can generate with xmaple-DIP

dot -Tpdf graph-$1.dot > uip.pdf; evince uip.pdf &
dot -Tpdf graph-current-$1.dot > uip2.pdf; evince uip2.pdf &
