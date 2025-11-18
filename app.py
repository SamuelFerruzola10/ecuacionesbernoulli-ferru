import sympy as sp
from flask import Flask, render_template, request, jsonify

x, C1 = sp.symbols('x C1') 

def solve_bernoulli_web(p_expr, q_expr, n_expr, ic=None):
    steps = []
    solucion_final = None
    solucion_particular = None
    p_latex = None
    q_latex = None

    try:
        p = sp.sympify(p_expr)
        q = sp.sympify(q_expr)
        n = sp.sympify(n_expr)
    except Exception as e:
        steps.append({"type": "error", "title": "Error de Entrada", "text": "Error al procesar las funciones p(x), q(x) o n: " + str(e)})
        return steps, None, None, None, None

    p_latex = sp.latex(p)
    q_latex = sp.latex(q)

    if sp.simplify(n) == 0 or sp.simplify(n - 1) == 0:
        n_val = "0" if sp.simplify(n) == 0 else "1"
        steps.append({"type": "error", "title": "Caso Trivial", "text": f"La ecuación es lineal ($n={n_val}$) y se resuelve con factor integrante, pero no aplica la sustitución de Bernoulli."})
        return steps, p_latex, q_latex, None, None

    steps.append({
        "type": "original", 
        "title": "Ecuación Original",
        "p": p_latex,
        "q": q_latex,
        "n": sp.latex(n),
        "formula": r"y' + p(x)y = q(x)y^n"
    })

    one_minus_n = sp.simplify(1 - n)
    Pu = sp.simplify(one_minus_n * p)
    R = sp.simplify(one_minus_n * q)
    
    steps.append({
        "type": "substitution", 
        "title": "Paso 1: Sustitución de Bernoulli", 
        "u": r"u = y^{" + sp.latex(one_minus_n) + "}",
        "P_u": sp.latex(Pu),
        "R": sp.latex(R),
        "formula": r"u' + P_u(x)u = R(x)"
    })

    try:
        IF = sp.exp(sp.integrate(Pu, x))
    except Exception as e:
        steps.append({"type": "error", "title": "Error Integrante", "text": "No se pudo calcular el Factor Integrante: " + str(e)})
        return steps, p_latex, q_latex, None, None

    steps.append({
        "type": "step", 
        "title": "Paso 2: Factor Integrante", 
        "text": r"IF = e^{\int P_u dx} = " + sp.latex(IF)
    })

    try:
        integral = sp.integrate(IF * R, x)
    except Exception as e:
        steps.append({"type": "error", "title": "Error Integral", "text": "No se pudo calcular la integral: " + str(e)})
        return steps, p_latex, q_latex, None, None

    steps.append({
        "type": "step", 
        "title": "Paso 3: Cálculo de la Integral", 
        "text": r"\int IF \cdot R dx = " + sp.latex(integral)
    })

    u_general = sp.simplify((integral + C1) / IF)
    steps.append({
        "type": "step", 
        "title": "Paso 4: Solución General para u(x)", 
        "text": r"u(x) = " + sp.latex(u_general)
    })

    exponent = sp.simplify(1 / one_minus_n)
    y_general_base = sp.simplify(u_general)
    y_general = r"\left(" + sp.latex(y_general_base) + r"\right)^{" + sp.latex(exponent) + r"}"
    
    steps.append({
        "type": "solution", 
        "title": "Paso 5: Solución General para y(x)", 
        "text": r"y(x) = " + y_general
    })
    
    solucion_final = y_general

    if ic is not None and len(ic) == 2:
        x0_val, y0_val = ic
        try:
            eq_for_C = sp.Eq(u_general.subs(x, x0_val), sp.sympify(y0_val) ** one_minus_n)
            solC = sp.solve(eq_for_C, C1)
            
            if solC:
                
                Csol_exacta = solC[0]
                Csol_redondeada = sp.N(Csol_exacta, 6) 
                
               
                u_particular = u_general.subs(C1, Csol_redondeada) 
                u_particular_simp = sp.simplify(u_particular)

                y_particular_base = sp.powsimp(u_particular_simp, force=True)
                solucion_particular = r"\left(" + sp.latex(y_particular_base) + r"\right)^{" + sp.latex(exponent) + r"}"

                steps.append({
                    "type": "step", 
                    "title": "Paso 6: Cálculo de la Constante", 
                    "text": f"Para $y({sp.latex(x0_val)}) = {sp.latex(y0_val)}$, la constante es $C_1 \\approx {sp.latex(Csol_redondeada)}$."
                })
                steps.append({
                    "type": "solution_particular", 
                    "title": "Paso 7: Solución Particular", 
                    "text": r"y(x) = " + solucion_particular
                })
        except Exception as e:
            steps.append({"type": "error_ci", "title": "Error CI", "text": "Error al calcular la constante: " + str(e)})

    return steps, p_latex, q_latex, solucion_final, solucion_particular


app = Flask(__name__)


@app.route('/')
def index():
    
    return render_template('ecua.html')

#@calculo
@app.route('/solve', methods=['POST'])
def solve():
    data = request.json
    p_in = data.get('p')
    q_in = data.get('q')
    n_in = data.get('n')
    x0_in = data.get('x0')
    y0_in = data.get('y0')
    
    ic = None
    if x0_in and y0_in:
        try:
            ic = (sp.sympify(x0_in), sp.sympify(y0_in))
        except Exception:
            pass 

    steps, p_latex, q_latex, general_sol_latex, particular_sol_latex = solve_bernoulli_web(p_in, q_in, n_in, ic=ic)
    
    return jsonify({
        'steps': steps,
        'p_latex': p_latex,
        'q_latex': q_latex,
        'n_value': n_in,
        'general_solution': general_sol_latex,
        'particular_solution': particular_sol_latex
    })

if __name__ == '__main__':
    app.run(debug=True)