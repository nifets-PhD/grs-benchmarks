import sys
import basico

basico.new_model(name="smoke-bd", volume=1.0, quantity_unit="#")
basico.add_species("X", initial_concentration=0)
basico.add_reaction("birth", "-> X")
basico.add_reaction("death", "X ->")
basico.set_reaction_parameters("(birth).v", value=10.0)
basico.set_reaction_parameters("(death).k1", value=0.1)

basico.set_task_settings(
    basico.T.TIME_COURSE,
    {
        "method": {
            "name": "Stochastic (Gibson + Bruck)",
            "Use Random Seed": True,
            "Random Seed": "1234",
        }
    },
)
res = basico.run_time_course(duration=200, stepnumber=5, method="stochastic")

print(
    f"copasi basico {basico.__version__}: birth death res = {res['X'].iloc[-1]:.0f}, expected ~100"
)
