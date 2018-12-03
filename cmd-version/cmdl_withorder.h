// ///////////////////////////// MIT License //////////////////////////////////// //
//                                                                                //
// Copyright (c) 2013 David Zsolt Manrique                                        //
//                    david.zsolt.manrique@gmail.com                              //
//                                                                                //
// Permission is hereby granted, free of charge, to any person obtaining a copy   //
// of this software and associated documentation files (the "Software"), to deal  //
// in the Software without restriction, including without limitation the rights   //
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell      //
// copies of the Software, and to permit persons to whom the Software is          //
// furnished to do so, subject to the following conditions:                       //
//                                                                                //
// The above copyright notice and this permission notice shall be included in     //
// all copies or substantial portions of the Software.                            //
//                                                                                //
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR     //
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,       //
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE    //
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER         //
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,  //
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN      //
// THE SOFTWARE.                                                                  //
//                                                                                //
// ////////////////////////////////////////////////////////////////////////////// //

#ifndef COMMAND_LINE_ARGUMENT_H
#define COMMAND_LINE_ARGUMENT_H

// standard headers
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <stdexcept>
#include <tuple>
#include <regex>

namespace cmdl {

    template<typename OFirstType>
    bool obtain(std::istream &is, OFirstType &ofirstvalue) {
        return is >> ofirstvalue;
    }

    template<typename OFirstType, typename... OTypes>
    bool obtain(std::istream &is, OFirstType &ofirstvalue, OTypes &... ovalues) {
        return (is >> ofirstvalue) && obtain(is, ovalues...);
    }

    class basearg {
        friend class arguments;

    protected:
        std::string _short_name;
        std::string _long_name;
        std::string _help_instruction;
        std::string _arg_type;
        int _input_number;
    };

//!
//!  command line argument class
//!  this contains the argc and argv, and the help message.
//!

    class arguments {
        template<typename T> friend
        class basetypearg;

    protected:
        std::string _command_line;
        std::string _program_name;
        std::list <std::string> _command_line_args;
        std::list <std::string> _original_command_line_args;

        std::vector<basearg *> _defined_args;

    public:

        arguments(int argc, char **argv) {
            std::ostringstream os;

            _program_name = argv[0];
            os << argv[0] << ' ';

            for (int i = 1; i < argc; i++) {
                _command_line_args.push_back(argv[i]);
                os << argv[i] << ' ';
            }
            _command_line_args.push_back("--");
            _command_line = os.str();
            _original_command_line_args = _command_line_args;
        }

        arguments(const std::string &cmdline) {
            _command_line = cmdline;

            std::istringstream iss(cmdline);
            std::vector <std::string> args;
            std::copy(std::istream_iterator<std::string>(iss),
                      std::istream_iterator<std::string>(),
                      std::back_inserter < std::vector < std::string > > (args));

            _program_name = args[0];
            for (int i = 1; i < args.size(); i++)
                _command_line_args.push_back(args[i]);
            _command_line_args.push_back("--");

            _original_command_line_args = _command_line_args;
        }

        const std::string &name() const {
            return _program_name;
        }

        const std::string &cmdline() const {
            return _command_line;
        }

        const std::list <std::string> &cmdline_args() const {
            return _command_line_args;
        }

        void print() const {
            for (std::list<std::string>::const_iterator it = _command_line_args.begin();
                 it != _command_line_args.end(); it++)
                std::cout << *it << std::endl;
            std::cout << std::endl;
        }

        //! print help messgae
        //! @param os stream

        void print_help(std::ostream &os) const {

            std::list <std::vector<std::string>> _help;
            for (std::vector<basearg *>::const_iterator it = _defined_args.begin(); it != _defined_args.end(); it++) {
                if ((*it)->_arg_type == "posarg") _help.push_back({(*it)->_help_instruction});
                else if ((*it)->_arg_type == "multiposarg") _help.push_back({(*it)->_help_instruction, "0"});
                else _help.push_back({(*it)->_short_name, (*it)->_long_name, (*it)->_help_instruction});
            }

            int n = 0;
            os << "Usage: " << _program_name << ' ';
            for (std::list < std::vector < std::string > > ::const_iterator it = _help.begin(); it != _help.end();
            it++)
            if (it->size() == 1) os << "<arg-" << ++n << ">" << ' ';
            else if (it->size() == 2 && (*it)[1] != "0") {
                os << "<arg-" << ++n << "...";
                std::istringstream is((*it)[1]);
                int dn = 0;
                is >> dn;
                n += dn;
                os << n << ">" << ' ';
            } else if (it->size() == 2 && (*it)[1] == "0") {
                os << "<arg-" << ++n << "... >" << ' ';
            }
            os << " -[option] <option-arg>" << std::endl;
            n = 0;
            for (std::list < std::vector < std::string > > ::const_iterator it = _help.begin(); it != _help.end();
            it++)
            if (it->size() == 1) {
                std::ostringstream oss;
                oss << "<arg-" << ++n << ">";
                os << "   " << std::left << std::setw(15) << oss.str() << ' ' << std::left << std::setw(50) << (*it)[0]
                   << std::endl;
            } else if (it->size() == 2 && (*it)[1] != "0") {
                std::ostringstream oss;
                oss << "<arg-" << ++n << "...";
                std::istringstream is((*it)[1]);
                int dn = 0;
                is >> dn;
                n += dn;
                oss << n << ">";
                os << "   " << std::left << std::setw(15) << oss.str() << ' ' << std::left << std::setw(50) << (*it)[0]
                   << std::endl;
            } else if (it->size() == 2 && (*it)[1] == "0") {
                std::ostringstream oss;
                oss << "<arg-" << ++n << "... >";
                os << "   " << std::left << std::setw(15) << oss.str() << ' ' << std::left << std::setw(50) << (*it)[0]
                   << std::endl;
            }
            for (std::list < std::vector < std::string > > ::const_iterator it = _help.begin(); it != _help.end();
            it++)
            if (it->size() == 3) {
                os << "   " << std::left << std::setw(5) << (*it)[0] << ' ';
                os << std::left << std::setw(20) << (*it)[1] << ' ';
                os << std::left << std::setw(50) << (*it)[2] << std::endl;
            }
        }

        void print_all() {
//            for (auto &da : _defined_args) {
//                std::cout << da->_arg_type << "\t" << da->_short_name << "\t" << da->_long_name << "\t" << da->_help_instruction << std::endl;
//            }
        }

        std::vector <std::pair<basearg *, size_t>> order() const {

            std::vector <std::pair<basearg *, size_t>> ord;
            std::map < basearg * , size_t > baseargmap;
            for (auto it = _original_command_line_args.begin(); it != _original_command_line_args.end(); ++it) {
                for (auto &da : _defined_args)
                    if (da->_short_name == *it || da->_long_name == *it) {
                        if (baseargmap.find(da) == baseargmap.end()) baseargmap[da] = 0;
                        ord.push_back(std::pair<basearg *, size_t>(da, baseargmap[da]));
                        baseargmap[da]++;
                    }
            }

            return ord;
        }
    };


    template<typename BaseType>
    class basetypearg : public basearg {
    protected:
        BaseType _var;
        bool _found;

        std::list <std::string> &_command_line_args;
        std::vector<basearg *> &_defined_args;
    public:

        basetypearg(arguments &args,
                    const std::string &short_name,
                    const std::string &long_name,
                    const std::string &help_instruction,
                    const std::string &arg_type, const int input_number) : _command_line_args(args._command_line_args),
                                                                           _defined_args(args._defined_args) {

            _defined_args.push_back(this);

            this->_short_name = short_name;
            this->_long_name = long_name;
            this->_help_instruction = help_instruction;
            this->_arg_type = arg_type;
            this->_input_number = input_number;

        }

        ~basetypearg() {
            _defined_args.erase(std::remove(_defined_args.begin(), _defined_args.end(), this), _defined_args.end());
        }

        const BaseType &operator*() const {
            return _var;
        }

        const BaseType &val() const {
            return _var;
        }

        bool is_set() const {
            return _found;
        }
    };

    template<typename Type>
    class posarg : public basetypearg<Type> {
    public:

        posarg(arguments &args, const std::string &help_instruction, Type initial_value) :
                basetypearg<Type>(args, "", "", help_instruction, "posarg", 1) {

            this->_var = initial_value;
            this->_found = false;

            for (std::list<std::string>::iterator it = this->_command_line_args.begin();
                 it != this->_command_line_args.end(); it++)
                if (*it != "--") {
                    std::istringstream is(*it);
                    if (!(is >> this->_var))
                        throw std::runtime_error("Wrong command line arguments. Try --help!");
                    it = this->_command_line_args.erase(it);
                    this->_found = true;
                    break;
                } else break;
        }
    };

    template<typename Type>
    class multiposarg : public basetypearg<std::vector < Type>

    > {
    public:

    multiposarg(arguments
    & args,
    const std::string &help_instruction, Type
    initial_value,
    const int length = 0
    ) :
    basetypearg<std::vector < Type> > (args, "", "", help_instruction, "multiposarg", length) {

    std::string str;
    std::stringstream ss;
    ss <<
    length;
    ss >>
    str;

    this->_var.

    clear();

    this->
    _found = false;

    for (
    std::list<std::string>::iterator it = this->_command_line_args.begin();
    it != this->_command_line_args.

    end();

    it++)
    if (*it != "--") {

    std::istringstream is(*it);

    Type tmp;

    if (!(is >> tmp))
    throw std::runtime_error("Wrong command line arguments. Try --help!");

    it = this->_command_line_args.erase(it);
    --
    it;

    this->_var.
    push_back(tmp);

    this->
    _found = true;

    if (length > 0 && this->_var.

    size()

    >= length) break;

}
else break;

if (!this->_found) this->_var.
push_back(initial_value);
}

};

template<bool default_value>
class switcharg : public basetypearg<bool> {
public:

    switcharg(arguments &args, const std::string &short_name, const std::string &long_name,
              const std::string &help_instruction) :
            basetypearg<bool>(args, short_name, long_name, help_instruction, "switcharg", 0) {

        this->_var = default_value;
        this->_found = false;

        for (std::list<std::string>::iterator it = this->_command_line_args.begin();
             it != this->_command_line_args.end(); it++)
            if (*it != "--") {
                if ((*it == short_name) || (*it == long_name)) {

                    if (this->_found)
                        throw std::runtime_error("Multiple " + *it + " command line arguments. Try --help!");

                    it = this->_command_line_args.erase(it);

                    this->_var = !_var;
                    this->_found = true;

                }
            } else break;
    }
};

template<typename Type>
class vararg : public basetypearg<Type> {
public:
    //! constructor
    //! @param args arguments type that contains the command line input
    //! @param s short version of the switch
    //! @param l long version of the switch
    //! @param help simple help message for this input
    //! @param def default value of this input

    vararg(arguments &args, const std::string &short_name, const std::string &long_name,
           const std::string &help_instruction, Type initial_value) :
            basetypearg<Type>(args, short_name, long_name, help_instruction, "vararg", 1) {

        this->_var = initial_value;
        this->_found = false;

        for (std::list<std::string>::iterator it = this->_command_line_args.begin();
             it != this->_command_line_args.end(); it++)
            if (*it != "--") {
                if ((*it == short_name) || (*it == long_name)) {

                    if (this->_found)
                        throw std::runtime_error("Multiple " + *it + " command line arguments. Try --help!");

                    if (*++it == "--")
                        throw std::runtime_error("Wrong command line arguments. Try --help!");

                    std::istringstream is(*it);
                    if (!(is >> this->_var))
                        throw std::runtime_error("Wrong command line arguments. Try --help!");

                    it = this->_command_line_args.erase(it);
                    it = this->_command_line_args.erase(--it);

                    this->_found = true;

                }
            } else break;

    }
};

template<typename FirstType, typename... Types>
class tuplevararg : public basetypearg<std::tuple < FirstType, Types...>

> {
public:

tuplevararg(arguments
& args,
const std::string &short_name,
const std::string &long_name,
const std::string &help_instruction, FirstType
firstvalue,
Types... values
) :
basetypearg <std::tuple<FirstType, Types...>> (args, short_name, long_name, help_instruction,
"tuplevararg", std::tuple_size<std::tuple < FirstType, Types...> >::value) {

this->
_var = std::tuple<FirstType, Types...>(firstvalue, values...);
this->
_found = false;

for (
std::list<std::string>::iterator it = this->_command_line_args.begin();
it != this->_command_line_args.

end();

it++)
if (*it != "--") {
if ((*it == short_name) || (*it == long_name)) {

if (this->_found)
throw std::runtime_error("Multiple " + *it + " command line arguments. Try --help!");

std::stringstream ios;

for (
int i = 0;
i<std::tuple_size < std::tuple < FirstType, Types...> >::value;
i++) {
if (*++it == "--")
throw std::runtime_error("Wrong command line arguments. Try --help!");
ios << *it << ' ';
}

if (!
obtain(ios, firstvalue, values
...))
throw std::runtime_error("Wrong command line arguments. Try --help!");

this->
_var = std::tuple<FirstType, Types...>(firstvalue, values...);

for (
int i = 0;
i<std::tuple_size < std::tuple < FirstType, Types...> >::value;
i++) {
it = this->_command_line_args.erase(it);
--
it;
}

it = this->_command_line_args.erase(it);
this->
_found = true;

}
} else break;
}
};

template<typename Type>
class multivararg : public basetypearg<std::vector < Type>

> {
public:

multivararg(arguments
& args,
const std::string &short_name,
const std::string &long_name,
const std::string &help_instruction, Type
initial_value) :
basetypearg <std::vector<Type>> (args, short_name, long_name, help_instruction,
"multivararg", 1) {

this->_var.

clear();

this->
_found = false;

for (
std::list<std::string>::iterator it = this->_command_line_args.begin();
it != this->_command_line_args.

end();

it++)
if (*it != "--") {
if ((*it == short_name) || (*it == long_name)) {
if (*++it == "--")
throw std::runtime_error("Wrong command line arguments. Try --help!");

std::istringstream is(*it);
Type tmp;

if (!(is >> tmp))
throw std::runtime_error("Wrong command line arguments. Try --help!");

this->_var.
push_back(tmp);

it = this->_command_line_args.erase(it);
it = this->_command_line_args.erase(--it);
it--;

this->
_found = true;
}
} else break;

if (!this->_found) this->_var.
push_back(initial_value);

}

};

template<typename FirstType, typename... Types>
class muplevararg : public basetypearg<std::vector < std::tuple < FirstType, Types...>

> > {
public:

muplevararg(arguments
& args,
const std::string &short_name,
const std::string &long_name,
const std::string &help_instruction, FirstType
firstvalue,
Types... values
) :
basetypearg <std::vector<std::tuple <
                         FirstType, Types...>> >(args, short_name, long_name, help_instruction, "muplevararg", std::tuple_size<
        std::tuple < FirstType, Types...> >::value) {

this->_var.

clear();

this->
_found = false;

for (
std::list<std::string>::iterator it = this->_command_line_args.begin();
it != this->_command_line_args.

end();

it++)
if (*it != "--") {
if ((*it == short_name) || (*it == long_name)) {

std::stringstream ios;

for (
int i = 0;
i<std::tuple_size < std::tuple < FirstType, Types...> >::value;
i++) {
if (*++it == "--")
throw std::runtime_error("Wrong command line arguments. Try --help!");
ios << *it << ' ';
}

if (!
obtain(ios, firstvalue, values
...))
throw std::runtime_error("Wrong command line arguments. Try --help!");

auto tmp = std::tuple<FirstType, Types...>(firstvalue, values...);

this->_var.
push_back(tmp);

for (
int i = 0;
i<std::tuple_size < std::tuple < FirstType, Types...> > ::value;
i++) {
it = this->_command_line_args.erase(it);
--
it;
}

it = this->_command_line_args.erase(it);
it--;

this->
_found = true;

}
} else break;

if (!this->_found) this->_var.
push_back(std::tuple<FirstType, Types...>(firstvalue, values...)
);
}

};

}
#endif // COMMAND_LINE_ARGUMENT_H

