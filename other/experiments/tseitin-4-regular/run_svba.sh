cp $1 ~/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin/s.cnf
cd ~/recerca/systems/SAT-comp-23-sequential/solvers/SBVA/bin && bash albert_run_sbva_cadical 40 20 s.cnf . > out.txt
