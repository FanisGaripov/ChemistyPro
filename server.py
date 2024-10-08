from flask import Flask, render_template, request, redirect, flash, session, url_for, jsonify, send_file
import re
from mod import db, User
from chempy import balance_stoichiometry
import os
import flask_login
import json
from flask_login import login_required, UserMixin, LoginManager, login_user
from werkzeug.utils import secure_filename
import requests
from bs4 import BeautifulSoup
from urllib.parse import quote
# импортируем все библиотеки

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///products.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = 'static/'
app.secret_key = 'supersecretkey'
login_manager = LoginManager()
login_manager.init_app(app)
db.init_app(app)
c = []

''' идеи:
1.Добавление таблицы растворимостей VVV
2.Ряд электроотрицательности элементов 
3.Электрохимический ряд активности металлов VV
4.Ряд силы кислот VV
5.Кислоты и кислотные остатки VV
'''

def molecular_mass(formula):
    # Словарь с атомными массами элементов
    atomic_masses = {
        'H': 1.008,
        'He': 4.0026,
        'Li': 6.94,
        'Be': 9.0122,
        'B': 10.81,
        'C': 12.011,
        'N': 14.007,
        'O': 15.999,
        'F': 18.998,
        'Ne': 20.180,
        'Na': 22.99,
        'Mg': 24.305,
        'Al': 26.982,
        'Si': 28.085,
        'P': 30.974,
        'S': 32.06,
        'Cl': 35.45,
        'Ar': 39.948,
        'K': 39.098,
        'Ca': 40.078,
        'Sc': 44.956,
        'Ti': 47.867,
        'V': 50.941,
        'Cr': 51.996,
        'Mn': 54.938,
        'Fe': 55.845,
        'Co': 58.933,
        'Ni': 58.693,
        'Cu': 63.546,
        'Zn': 65.38,
        'Ga': 69.723,
        'Ge': 72.630,
        'As': 74.922,
        'Se': 78.971,
        'Br': 79.904,
        'Kr': 83.798,
        'Rb': 85.468,
        'Sr': 87.62,
        'Y': 88.906,
        'Zr': 91.224,
        'Nb': 92.906,
        'Mo': 95.95,
        'Tc': 98,
        'Ru': 101.07,
        'Rh': 102.905,
        'Pd': 106.42,
        'Ag': 107.868,
        'Cd': 112.414,
        'In': 114.818,
        'Sn': 118.710,
        'Sb': 121.760,
        'Te': 127.60,
        'I': 126.904,
        'Xe': 131.293,
        'Cs': 132.905,
        'Ba': 137.327,
        'La': 138.905,
        'Ce': 140.116,
        'Pr': 140.907,
        'Nd': 144.242,
        'Pm': 145,
        'Sm': 150.36,
        'Eu': 151.964,
        'Gd': 157.25,
        'Tb': 158.925,
        'Dy': 162.500,
        'Ho': 164.930,
        'Er': 167.259,
        'Tm': 168.934,
        'Yb': 173.04,
        'Lu': 174.966,
        'Hf': 178.49,
        'Ta': 180.947,
        'W': 183.84,
        'Re': 186.207,
        'Os': 190.23,
        'Ir': 192.217,
        'Pt': 195.084,
        'Au': 196.967,
        'Hg': 200.592,
        'Tl': 204.38,
        'Pb': 207.2,
        'Bi': 208.980,
        'Po': 209,
        'At': 210,
        'Rn': 222,
        'Fr': 223,
        'Ra': 226,
        'Ac': 227,
        'Th': 232.038,
        'Pa': 231.035,
        'U': 238.028,
        'Np': 237,
        'Pu': 244,
        'Am': 243,
        'Cm': 247,
        'Bk': 247,
        'Cf': 251,
        'Es': 252,
        'Fm': 257,
        'Md': 258,
        'No': 259,
        'Lr': 262,
        'Rf': 267,
        'Db': 270,
        'Sg': 271,
        'Bh': 270,
        'Hs': 277,
        'Mt': 276,
        'Ds': 281,
        'Rg': 282,
        'Cn': 285,
        'Nh': 286,
        'Fl': 289,
        'Mc': 290,
        'Lv': 293,
        'Ts': 294,
        'Og': 294,
    }

    mass = 0.0
    elements = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    element_details = []  # Список для хранения деталей элементов
    for element, count in elements:
        count = int(count) if count else 1  # Установите 1, если count отсутствует
        element_mass = atomic_masses[element]
        total_mass = element_mass * count
        element_details.append((element, atomic_masses[element], count, total_mass))  # Добавляем элемент, его массу, количество и общую массу
        mass += total_mass

    return round(mass), element_details  # Возвращаем общую массу и детали элементов


def electronic_configuration(element):
    elements_data = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19,
        'Ca': 20, 'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28,
        'Cu': 29, 'Zn': 30, 'Ga': 31, 'Ge': 32, 'As': 33, 'Se': 34, 'Br': 35, 'Kr': 36, 'Rb': 37,
        'Sr': 38, 'Y': 39, 'Zr': 40, 'Nb': 41, 'Mo': 42, 'Tc': 43, 'Ru': 44, 'Rh': 45, 'Pd': 46,
        'Ag': 47, 'Cd': 48, 'In': 49, 'Sn': 50, 'Sb': 51, 'Te': 52, 'I': 53, 'Xe': 54, 'Cs': 55,
        'Ba': 56, 'La': 57, 'Ce': 58, 'Pr': 59, 'Nd': 60, 'Pm': 61, 'Sm': 62, 'Eu': 63, 'Gd': 64,
        'Tb': 65, 'Dy': 66, 'Ho': 67, 'Er': 68, 'Tm': 69, 'Yb': 70, 'Lu': 71, 'Hf': 72, 'Ta': 73,
        'W': 74, 'Re': 75, 'Os': 76, 'Ir': 77, 'Pt': 78, 'Au': 79, 'Hg': 80, 'Tl': 81, 'Pb': 82,
        'Bi': 83, 'Po': 84, 'At': 85, 'Rn': 86, 'Fr': 87, 'Ra': 88, 'Ac': 89, 'Th': 90, 'Pa': 91,
        'U': 92, 'Np': 93, 'Pu': 94, 'Am': 95, 'Cm': 96, 'Bk': 97, 'Cf': 98, 'Es': 99, 'Fm': 100,
        'Md': 101, 'No': 102, 'Lr': 103, 'Rf': 104, 'Db': 105, 'Sg': 106, 'Bh': 107, 'Hs': 108,
        'Mt': 109, 'Ds': 110, 'Rg': 111, 'Cn': 112, 'Nh': 113, 'Fl': 114, 'Mc': 115, 'Lv': 116,
        'Ts': 117, 'Og': 118
    }

    atomic_number = elements_data.get(element)
    if element == '':
        return "Введите элемент"
    if atomic_number is None:
        return "Элемент не найден"

    configurations = []
    '''subshells = ['1s', '2s', '2p', '3s', '3p', '4s', '3d', '4p', '5s', '4d', '5p', '6s', '4f', '5d', '6p', '7s', '5f',
                 '6d', '7p']
    electrons = [2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6]'''
    subshells = ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '4f', '5s', '5p', '5d', '5f', '6s', '6p', '6d',
                 '7s', '7p']
    electrons = [2, 2, 6, 2, 6, 10, 2, 6, 10, 14, 2, 6, 10, 14, 2, 6, 10, 2, 6]

    for i in range(len(subshells)):
        if atomic_number > 0:
            if atomic_number >= electrons[i]:
                configurations.append(f"{subshells[i]}^{electrons[i]}")
                atomic_number -= electrons[i]
            else:
                configurations.append(f"{subshells[i]}^{atomic_number}")
                break
    return ' '.join(configurations)


@app.route('/electronic_configuration', methods=['GET', 'POST'])
def electronic_configuration_page():
    # функция, которая отображает страницу электронной конфигурации, предыдущая функция отвечает за обработку ответа
    element = ''
    user = flask_login.current_user
    configuration = ''
    if request.method == 'POST':
        element = request.form.get("element", False)
        configuration = electronic_configuration(element)
    return render_template('electronic_configuration.html', configuration=configuration, user=user, element=element)


def uravnivanie(formula):
    # баланс уравнений
    reactants_input, products_input = formula.split('=')
    reactants = {x.split()[0].strip(): int(x.split()[1]) if len(x.split()) > 1 else 1 for x in
                 reactants_input.split('+')}
    products = {x.split()[0].strip(): int(x.split()[1]) if len(x.split()) > 1 else 1 for x in products_input.split('+')}

    balanced_reaction = balance_stoichiometry(reactants, products)

    reactants_str = ' + '.join([f"{v}{k}" for k, v in balanced_reaction[0].items()])
    products_str = ' + '.join([f"{v}{k}" for k, v in balanced_reaction[1].items()])

    otvet = f"{reactants_str} = {products_str}"

    return otvet


@app.route('/', methods=['GET', 'POST'])
def osnova():
    # функция которая возвращает главную страницу сайта( index.html )
    user = flask_login.current_user
    resultat2 = ''
    if request.method == 'POST':
        chemical_formula = request.form['chemical_formula']
        try:
            resultat2 = f'{chemical_formula}: {uravnivanie(chemical_formula)}'
            print(resultat2)
        except:
            redirect('/')
    return render_template('index.html', resultat2=resultat2, user=user)


@app.route('/molyarnaya_massa', methods=['GET', 'POST'])
def molyar_massa():
    # метод для вычисления молярной массы и отображения ее на сайте
    user = flask_login.current_user
    global resultat, dlyproverki, c
    resultat = ''
    otdelno = []
    formatspisok = ''
    dlyproverki = 0
    if request.method == 'POST':
        chemical_formula = request.form['element']
        try:
            dlyproverki, element_details = molecular_mass(chemical_formula)
            resultat = f"Молярная масса {chemical_formula}: {dlyproverki} г/моль"
            for element, mass, count, total_mass in element_details:
                otdelno.append(f"{count} x {element} ({round(mass)} г/моль): {round(total_mass)} г")
        except Exception as e:
            print(f"Ошибка: {e}")
            return redirect('/')
    return render_template('molyarnaya_massa.html', resultat=resultat, dlyproverki=dlyproverki, user=user, otdelno=otdelno)


def get_chemical_equation_solution(reaction):
    if request.method == 'POST':
        reaction = request.form.get("chemical_formula", False)
    # Кодируем реакцию для URL
        encoded_reaction = quote(reaction)

    # Формируем URL с учетом химической реакции
        url = f"https://chemequations.com/ru/?s={encoded_reaction}"

    # Отправляем GET-запрос
        response = requests.get(url)

    # Проверка успешности запроса
        if response.status_code == 200:
            # Парсим HTML-ответ
            soup = BeautifulSoup(response.text, 'html.parser')

            # Находим элемент с классом "equation main-equation well"
            result = soup.find('h1', class_='equation main-equation well')

            if result:
                return result.get_text(strip=True)
                # Возвращаем текст ответа
            else:
                return 'Решение не найдено.'
        else:
            return f"Ошибка при запросе: {response.status_code}"


@app.route('/complete_reaction', methods=['GET', 'POST'])
def complete_reaction_page():
    # страница, отвечающая за вывод завершенных реакций предыдущим методом
    react1 = ''
    user = flask_login.current_user
    reaction = ''
    if request.method == 'POST':
        reaction = request.form.get("chemical_formula", False)
        react1 = get_chemical_equation_solution(reaction)

    return render_template('complete_reaction.html', get_chemical_equation_solution=get_chemical_equation_solution, react1=react1, user=user, reaction=reaction)


def get_reaction_chain(reaction):
    # цепочка превращений
    if request.method == 'POST':
        reaction = request.form.get("chemical_formula", False)
        url = f"https://chemer.ru/services/reactions/chains/{reaction}"
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3'
        }
        session = requests.Session()
        session.headers.update(headers)
        response = session.get(url)
        if response.status_code == 200:
            soup = BeautifulSoup(response.text, 'html.parser')
            content_sections = soup.find_all('section', class_='content')  # Ищем все секции с классом 'content'

            results = []  # Список для хранения первых ответов из каждой секции

            for content_section in content_sections:
                reactions = content_section.find_all('p', class_='resizable-block')  # Ищем все 'p' внутри каждой секции
                if reactions:
                    first_reaction = reactions[0].get_text().strip()  # Берем первый ответ из секции
                    results.append(first_reaction)  # Добавляем его в результаты

            return results  # Возвращаем список первых ответов из всех секций


@app.route('/get_reaction_chain', methods=['GET', 'POST'])
def get_reaction_chain_page():
    # страница, которая выводит цепочку превращений, т.е прошлую функцию
    user = flask_login.current_user
    react2 = ''
    reaction = ''
    if request.method == 'POST':
        reaction = request.form.get("chemical_formula", False)
        react2 = get_reaction_chain(reaction)

    return render_template('get_reaction_chain.html', get_reaction_chain=get_reaction_chain, user=user, reaction=reaction, react2=react2)


@app.route('/aboutme', methods=['GET', 'POST'])
def aboutme():
    # обо мне
    user = flask_login.current_user
    if user.is_authenticated:
        return render_template('about.html', user=user)
    else:
        return render_template('login.html', user=user)


@app.route('/instruction', methods=['GET', 'POST'])
def instruction():
    # инструкция
    user = flask_login.current_user
    return render_template('instruction.html', user=user)


@app.route('/tablica', methods=['GET', 'POST'])
def tablica():
    # таблица менделеева
    user = flask_login.current_user
    return render_template('tablica.html', user=user)


@app.route('/tablica_rastvorimosti', methods=['GET', 'POST'])
def tablica_rastvorimosti():
    # таблица растворимости
    user = flask_login.current_user
    return render_template('tablica_rastvorimosti.html', user=user)


@app.route('/tablica_kislotnosti', methods=['GET', 'POST'])
def tablica_kislotnosti():
    # таблица кислот ( ошибка в названии функции ) :))
    user = flask_login.current_user
    return render_template('tablica_kislotnosti.html', user=user)


@app.route('/login', methods=['GET', 'POST'])
def login():
    # вход пользователя в аккаунт. берутся данные из базы данных
    user = flask_login.current_user
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        user = User.query.filter_by(username=username).first()
        if user and user.check_password(password):
            login_user(user)
            return redirect('/')
        else:
            return render_template('bug.html', user=user)
    return render_template('login.html', user=user)


@login_manager.user_loader
def load_user(user_id):
    return User.query.get(int(user_id))


@login_manager.request_loader
def load_user_from_request(request):
    user_id = request.args.get('user_id')
    if user_id:
        return User.query.get(user_id)
    return None


@app.route('/register', methods=['GET', 'POST'])
def register():
    # регистрация, идет работа с бд
    user = flask_login.current_user
    if request.method == 'POST':
        username = request.form['username']
        password = request.form['password']
        name = request.form['name']
        surname = request.form['surname']
        email = request.form['email']
        if not User.query.filter_by(username=username).first():
            user = User(username=username, name=name, surname=surname, email=email)
            user.set_password(password)
            db.session.add(user)
            db.session.commit()
            return redirect(url_for('login'))
        else:
            flash('Пользователь с таким именем уже существует.')
    return render_template('register.html', user=user)


@app.route('/profile')
def profile():
    # профиль
    user = flask_login.current_user
    if user.is_authenticated:
        return render_template('profile.html', user=user)
    else:
        return redirect(url_for('login'))


with app.app_context():
    db.create_all()


@app.route('/logout')
@login_required
def logout():
    flask_login.logout_user()
    return redirect(url_for('login'))


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
