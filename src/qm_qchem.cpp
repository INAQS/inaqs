#include <ctime>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdexcept>
#include <sys/stat.h>
#include <algorithm>
#include "qm_qchem.hpp"
#include "properties.hpp"
#include "util.hpp"
#include <armadillo>
#include <unordered_map>
#include <chrono>
#include <functional>

ConfigBlockReader
QM_QChem::qchem_reader() {
  using types = ConfigBlockReader::types;
  ConfigBlockReader reader{"qchem"};
  reader.add_entry("qc_scratch", "DEFAULT");
  reader.add_entry("basis", types::STRING);
  reader.add_entry("verbatim_file", ""); // sentinel value
  reader.add_entry("exchange", types::STRING);
  reader.add_entry("scf_algorithm", "DIIS");
  reader.add_entry("scf_convergence", 8);
  reader.add_entry("two_electron_thresh", 11);
  reader.add_entry("scf_guess", "read");
  reader.add_entry("diabatization_method", "boys");
  reader.add_entry("diabat_states", std::vector<int> {}); // sentinel value
  reader.add_entry("donor_acceptor_ref", "");  // sentinel value
  reader.add_entry("nthreads", 1);
  reader.add_entry("buffer_states", 0);
  reader.add_entry("spin_diabats", std::vector<int> {});  // sentinel value

  reader.add_entry("sing_thresh", 1.2); // for state tracking

  // FIXME: ConfigReader should support bool
  reader.add_entry("singlets", 1);
  reader.add_entry("triplets", 0);
  reader.add_entry("unrestricted", 1);
  reader.add_entry("state_analysis", 0);
  reader.add_entry("spin_flip", 0);
  reader.add_entry("record_spectrum", 0);
  reader.add_entry("track_states", 0);
  reader.add_entry("strictly_diabatic_approximation", 1);
  reader.add_entry("loc_cis_ov_separate", 0);
  reader.add_entry("dump_qc_output", 0);
  reader.add_entry("dump_qc_input", 0);
  reader.add_entry("externalcharges_hack", 0);
  reader.add_entry("time_properties", 0);

  return reader;
}


/*
  Recall that Q-Chem can decide to change the number of excited states
  it computes! Usually this will be to a larger number of states. --Y.S.

  Happens in setman; to see where, search for "NRoots was altered as:"

  On subsequent invocation (as during AIMD), setman_init.F will reset
  the number requested to the original number.

  This is fixed by setting set_roots_orig and making sure we have
  headroom on our hopping states.
*/
QM_QChem::QM_QChem(FileHandle& fh, 
                   const arma::uvec& in_qmids, 
		   arma::mat& in_qm_crd, 
		   arma::mat& in_mm_crd, 
		   arma::vec& in_mm_chg, 
		   const int charge, 
		   const int mult,
		   const int excited_states_in,
                   const int min_state_in):
  QMInterface(in_qmids, in_qm_crd, in_mm_crd, in_mm_chg, charge, mult, excited_states_in, min_state_in)
{
  ConfigBlockReader reader = qchem_reader();
  reader.parse(fh);

  reader.get_data("basis", basis_set);
  reader.get_data("verbatim_file", verbatim_file);
  reader.get_data("exchange", exchange_method);
  reader.get_data("scf_algorithm", scf_algorithm);
  reader.get_data("scf_convergence", scf_convergence);
  reader.get_data("two_electron_thresh", two_electron_thresh);
  reader.get_data("scf_guess", scf_guess);
  reader.get_data("diabatization_method", diabatization_method);
  reader.get_data("sing_thresh", sing_thresh);

  // FIXME: maybe do this maybe not?
  // if (two_electron_thresh < scf_convergence + 3){
  //   std::cerr << "[QM_QChem] " << call_idx() << ": WARNING: two_electron_thresh should generally be 3 higher than scf_convergece." << std::endl;
  // }
  
  {
    int in;
    reader.get_data("singlets", in);
    singlets = in;
    reader.get_data("triplets", in);
    triplets = in;
    reader.get_data("unrestricted", in);
    unrestricted = in;    
    reader.get_data("state_analysis", in);
    state_analysis = in;
    reader.get_data("spin_flip", in);
    spin_flip = in;
    reader.get_data("record_spectrum", in);
    record_spectrum = in;
    reader.get_data("track_states", in);
    track_states = in;
    reader.get_data("strictly_diabatic_approximation", in);
    strictly_diabatic_approximation = in;
    reader.get_data("loc_cis_ov_separate", in);
    loc_cis_ov_separate = in;
    reader.get_data("dump_qc_output", in);
    dump_qc_output = in;
    reader.get_data("dump_qc_input", in);
    dump_qc_input = in;
    reader.get_data("externalcharges_hack", in);
    externalcharges_hack = in;
    reader.get_data("time_properties", in);
    time_properties = in;

  }

  // shstates must be set before buffer states or are added or min_state is changed for spin_flip
  shstates = excited_states + 1 - min_state;

  if(spin_flip){
    qm_multiplicity = 3;
    min_state += 1; // since the ground state is the first excited state
    excited_states += 1;
  }

  bool valid_options = true;

  if (min_state < 1){
    std::cerr << "min_state=" << min_state << std::endl;
    std::cerr << "Minimum states smaller than 1 not supported in QChem!" << std::endl;
    valid_options = false;
  }

  {
    int buffer_states;
    reader.get_data("buffer_states", buffer_states);

    if (buffer_states < 0){
      valid_options = false;
      std::cerr << "Buffer_states must be non-negative" << std::endl;
    }
    else{
      excited_states += buffer_states;
    }
  }

  {
    std::string in;
    reader.get_data("donor_acceptor_ref", in);
    if (in.length() > 0){
      arma::mat out(in);
      if (out.n_cols == 3 && out.n_rows == 2){
        donor_acceptor_ref = out.t();
        donor_acceptor_ref.print("Reference DA dipoles:");
      }
      else{
        valid_options = false;
        std::cerr << "Invalid attempt to set donor_acceptor_ref; must be of the form:" << std::endl;
        std::cerr << "  donor_acceptor_ref = Dx Dy Dz ; Ax Ay Az "<< std::endl;
        std::cerr << "where Di and Ai are respectively the dipole moments for the donor and acceptor states." << std::endl;
      }
    }
  }

  //FIXME: need to do better error checking
  reader.get_data("spin_diabats", spin_diabats);
  for (auto m: spin_diabats){
      if (m < 1){
        valid_options = false;
        std::cerr << "Invalid multiplicity for spin diabats: m=" << m << "!" << std::endl;
      }
    }  
  
  reader.get_data("diabat_states", diabat_states);
  if (diabat_states.size() > 1){
    //FIXME: when we have capacity for diabatic dynamics, will need to expand control
    boys_diabatization = true;

    for (auto s: diabat_states){
      if ((s < 0 ) || ((size_t) s > excited_states)){
        valid_options = false;
        std::cerr << "Diabatization states must be within the excited states!" << std::endl;
      }
    }

    // rely on `&&` to short-circut
    if (spin_diabats.size() > 0 && spin_diabats.size() != diabat_states.size()){
      valid_options = false;
      std::cerr << "If specifying both diabat_states and spin_diabats, they must be of the same lenght!" << std::endl;
    }
  }

  if(spin_flip && excited_states < 1){
    valid_options = false;
    std::cerr << "Must compute excited states with spin-flip; try increasing buffer_states." << std::endl;
  }

  if(spin_flip && !track_states){
    //valid_options = false;
    std::cerr << "WARNING: Selecting spin-flip without state tracking nearly guaranteed to produce incorrect results" << std::endl;
  }

  if (singlets && triplets && !track_states){
    valid_options = false;
    std::cerr << "Selecting singlets and triplets without state tracking nearly guaranteed to produce incorrect results!" << std::endl;
  }

  // FIXME: implement triplet tracking
  if (track_states && !singlets){
    valid_options = false;
    std::cerr << "State tracking for triplets not (yet << std::endl implemented!" << std::endl;
  }

  if (!valid_options){
    throw std::runtime_error("Invalid options; see above!");
  }
  
  /*
    The environmental variables $QCSCRATCH, $QCTHREADS shall take
    precedence over anything set in the config.
  */

  std::string conf_scratch;
  reader.get_data("qc_scratch", conf_scratch);
  qc_scratch_directory = get_qcscratch(conf_scratch);
  
  /*
    Setting $QCTHREADS seems sufficient on the Subotnik cluster, but
    perhaps this is unique to our setup. Need to check with Evgeny to
    be sure. Changes to $OMP_NUM_THREADS seem to have no effect.
  */
  int conf_threads;
  reader.get_data("nthreads", conf_threads);
  if (conf_threads > 1){
    setenv("QCTHREADS", std::to_string(conf_threads).c_str(), 0);
  }
  // FIXME: should include MPI support

  qc_executable = get_qcprog();
}


/*
  FIXME:
    - unify qm & mm gradients (all properties) [okay for link atoms?]
    - unify *_multi so we just give the indexes and a cube if we want multiple [why not use our kludge polymorphism]
*/
void QM_QChem::get_properties(PropMap &props){
  if(track_states){
    state_tracker(props);
  }

  for (QMProperty p: props.keys()){
    auto time_start = std::chrono::steady_clock::now();

    switch(p){
    case QMProperty::wfoverlap:{
      get_wf_overlap(props.get(p));
      break;
    }

    // FIXME: The scheme for diabats is horribly coupled
    case QMProperty::diabatic_rot_mat: break; // if we need U or H, pick them up with Ga
    case QMProperty::diabatic_H: break;
    case QMProperty::diabatic_gradients:{
      size_t A = 0, B = 0;

      const arma::uvec &idx = *props.get_idx(p);
      if (idx.size() > 1){
        A = idx[0];
        B = idx[1];
      }
      else if (diabat_states.size() > 1){
        A = diabat_states[0];
        B = diabat_states[1];
      }
      

      arma::mat H(2, 2, arma::fill::zeros), U(2, 2, arma::fill::zeros);
      arma::cube & gd = props.get(p);
      // FIXME: should set size for ecerything that asks for it!
      gd.set_size(3, NQM, 2);
      get_diabats(gd, U, H, A, B);

      if (props.has(QMProperty::diabatic_rot_mat)){
        props.get(QMProperty::diabatic_rot_mat) = U;
      }
      if (props.has(QMProperty::diabatic_H)){
        props.get(QMProperty::diabatic_H) = H;
      }
      
      break;
    }
      
    case QMProperty::nacvector_imag:
      throw std::invalid_argument("Imaginary NAC not implemented!");
      break;

    case QMProperty::nacvector:{
      const arma::uvec &idx = *props.get_idx(p);
      size_t I = idx[0]; size_t J = idx[1];
      get_nac_vector(props.get(p), I, J);
      break;
    }

    case QMProperty::mmgradient:
      break; // if we need mm, then we will do qm, which catches both
    case QMProperty::qmgradient:{
      arma::mat &g_qm = props.get(QMProperty::qmgradient);

      arma::uword surface = 0; //assume ground state
      if (props.has_idx(QMProperty::qmgradient)){
	surface = (*props.get_idx(QMProperty::qmgradient))[0];
      }

      if (props.has(QMProperty::mmgradient)){
	arma::mat &g_mm = props.get(QMProperty::mmgradient);
	get_gradient(g_qm, g_mm, surface);
      }
      else{
	get_gradient(g_qm, surface);
      }
      break;
    }

    case QMProperty::mmgradient_multi:
      break; // if we need mm, then we will do qm, which catches both
    case QMProperty::qmgradient_multi:{
      arma::cube &g_qm = props.get(QMProperty::qmgradient_multi);

      if (!props.has_idx(QMProperty::qmgradient_multi)){
        throw std::logic_error("Cannot request gradient_multi without specifying states!");
      }
      arma::uvec surfaces = *props.get_idx(QMProperty::qmgradient_multi);

      if (g_qm.n_slices != surfaces.n_elem){
	throw std::range_error("Insufficient space for requested gradients!");
      }
      
      for(arma::uword i = 0; i < surfaces.n_elem; i++){
	if (props.has(QMProperty::mmgradient_multi)){
	  arma::cube &g_mm = props.get(QMProperty::mmgradient_multi);
	  get_gradient(g_qm.slice(i), g_mm.slice(i), surfaces(i));
	}
	else{
	  get_gradient(g_qm.slice(i), surfaces(i));
	}
      }
      break;
    }
      
    case QMProperty::energies:{
      arma::vec & energies = props.get(p);
      if (props.has_idx(p)){
        arma::uvec idx = *props.get_idx(p);
        if (energies.n_elem != idx.n_elem){
          throw std::logic_error("energy array of incorrect size");
        }

        arma::vec e_temp(excited_states + 1, arma::fill::zeros);
        get_all_energies(e_temp);
        energies = e_temp(idx);
      }
      else{
        get_ground_energy(energies); // FIXME: combine ground/excited calls
      }
      break;
    }

    default: // Compile with the default case commented out and the compiler will detect neglected properties
      std::cerr << "You requested QMProperty " << p << std::endl;
      throw std::invalid_argument("Unknown QMProperty!");
      break;
    }
    auto time_stop = std::chrono::steady_clock::now();
    if (time_properties){
      std::cout << p << " required " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(time_stop - time_start).count()/1e3 << "s" << std::endl;
    }
  }

  // everything here will be called exactly once per call to update()
  if (!called(Q::once)){
    auto time_start = std::chrono::steady_clock::now();
    if (state_analysis){
      do_state_analysis();
    }

    if (record_spectrum){
      do_record_spectrum();
    }

    if (boys_diabatization){
      do_boys_diabatization();
    }

    auto time_stop = std::chrono::steady_clock::now();
    if (time_properties){
      std::cout << "once calls required " <<
        std::chrono::duration_cast<std::chrono::milliseconds>(time_stop - time_start).count() / 1e3 << "s" << std::endl;
    }
  }
}


// modifies all properties requesting state k to request k' where k' is the (k+1)th *singlet*
// Consider hooking into QC's state tracker via FILE_MECP_INFO
void QM_QChem::state_tracker(PropMap &props){
  arma::vec S2(excited_states + 1);
  {
    get_all_energies(S2); // a dummy call to generate FILE_CIS_S2; incidentally will also cache results for us
    S2.set_size(excited_states);
    readQFMan(QCFILE::FILE_CIS_S2, S2);
  }

  S2.save(arma::hdf5_name("gifs.hdf5",
                          "/spin2/" + std::to_string(call_idx()),
                          arma::hdf5_opts::replace));


  /* 1.2 is the default for REM_CIS_S2_THRESH; If we want to change,
     we can specify our threshold via sing_thresh, which is synced
     with the rem section.
  */

  // if S^2 is larger than sing_thresh, we consider the state a triplet

  arma::uvec statei(excited_states + 1, arma::fill::zeros);

  // FIXME: once BOMD doesn't know anything about excited states or
  // min_state (and idx{0} means the first comptued singlet, then we
  // can ditch these hacks. rather could be: n_singlets = 0;
  // statei(n_singlets++) = i; And we will need to test the ground
  // state too.
  arma::uword n_singlets = qm_multiplicity == 1 ? 1 : 0;

  for (arma::uword i = 0; i < excited_states; i++){
    if (S2(i) < sing_thresh){
      statei(n_singlets++) = i+1; // the ground state is state 0
    }
  }

  if (n_singlets == 1 && excited_states > 0){
    throw std::runtime_error("No states of the desired multiplicity computed!");
  }
  
  statei.resize(n_singlets);

  for (QMProperty p: props.keys()){
    if (props.has_idx(p)){
      arma::uvec & idxs = props.get_writable_idx(p);
      for (arma::uword & i: idxs){
        if (i >= n_singlets){
          std::cerr << "Requested " << p << "(" << i << "); have up to" << n_singlets-1 << "." << std::endl;
          throw std::runtime_error("The requested state is not within those states computed!");
        }
        else{
          i = statei(i);
        }
      }
    }
  }
}


void QM_QChem::do_boys_diabatization(void){
  REMKeys keys = {{"jobtype","sp"},
                  {diabatization_method + "_cis_numstate", std::to_string(diabat_states.size())},
                  {"loc_cis_ov_separate", std::to_string(loc_cis_ov_separate)},
                  {"sts_mom", "1"}};

  REMKeys ex = excited_rem();
  keys.insert(ex.begin(), ex.end());

  std::ofstream input = get_input_handle();
  write_rem_section(input, keys);  
  write_molecule_section(input);

  { // write boys section
    input << "$localized_diabatization" << std::endl;
    input << "Comment: states we mix" << std::endl;
    for (const auto& s: diabat_states){
      input << s << " ";
    }
    input << std::endl << "$end" << std::endl;
  }
  input.close();

  exec_qchem();

  // FIXME: write the relevant boys quantities in a better fashion
  {
    std::string src = get_qcwdir() + "/" + qc_log_file;
    std::string dest = get_qcwdir() + "/" +
      "boys." + std::to_string(call_idx()) + ".out";

    std::string cmd = "cp -a " + src + " " + dest; 
    int status = std::system(cmd.c_str());

    if(status){
      std::cerr << "Warning: unable to write " + dest << std::endl;
    }
  }
}

std::string QM_QChem::qc_log_file_idx(void){
  const int idx = call_idx();
  static int last_idx = idx;
  static int call = 0; // will be incremented on first run

  if (idx == last_idx){
    call++;
  }
  else{
    call = 1;
  }

  last_idx = idx;
  return "QCHEM." + std::to_string(idx) + "." + std::to_string(call) + ".out";
}


void QM_QChem::do_record_spectrum(void){
  arma::vec e;
  // dummy call to make sure the transition dipoles and energies are written
  e.set_size(excited_states + 1);
  get_all_energies(e);

  e.set_size(excited_states);
  readQFMan(QCFILE::FILE_SET_ENERGY, e);
  
  // for mu, the first row is the oscillator strength
  arma::mat mu(4, excited_states);
  readQFMan(QCFILE::FILE_TRANS_DIP_MOM, mu);

  //write spectrum to file as [excitation energy] [strength]
  {
    if (e[0] < 0){ // the true ground state is the first excited state
      e = e - e[0];
    }
    
    arma::mat spec(excited_states,2);
    spec.col(0) = e;
    spec.col(1) = mu.row(0).t();

    {
      arma::vec osc = spec.col(1);
      osc.save(arma::hdf5_name("gifs.hdf5",
                               "/oscillator/" + std::to_string(call_idx()),
                               arma::hdf5_opts::replace));
    }

    {
      std::string specf = get_qcwdir() + "/" + "spectrum.dat";
      std::ofstream stream;
      stream.open(specf, std::ios::out | std::ios::app | std::ios::binary);
      if (!stream){
        throw std::runtime_error("Cannot open spectrum file, " + specf);
      }
      else{
        spec.save(stream, arma::raw_ascii);
      }
    }
  }

}


void QM_QChem::do_state_analysis(void){
  REMKeys keys = {{"jobtype","sp"},
                  {"make_cube_files", "false"},
                  {"gui", "2"},
                  {"state_analysis", "true"},
                  {"molden_format", "true"}};

  REMKeys ex = excited_rem();
  keys.insert(ex.begin(), ex.end());

  std::ofstream input = get_input_handle();
  write_rem_section(input, keys);  
  write_molecule_section(input);
  input.close();

  {
    std::string guifile = get_qcwdir() + "/" +
      std::to_string(call_idx()) + ".fchk";
    setenv("GUIFILE", guifile.c_str(), 1);
  }
  
  exec_qchem();
}


//FIXME: Use another REM variable to save PREV_GEO
void QM_QChem::get_wf_overlap(arma::mat &U){
  static bool first_wfoverlap = true;

  if (U.n_elem != shstates * shstates){
    throw std::logic_error("WF overlap is of wrong size!");
  }
  
  REMKeys keys = excited_rem();
  keys.insert({{"jobtype","sp"},
      //{"namd_lowestsurface", std::to_string(min_state)},
               {"wf_overlap_minsurf", std::to_string(min_state)},
               {"wf_overlap_nsurf", std::to_string(shstates)},
               {"dump_wf_overlap", std::to_string(first_wfoverlap ? 1 : 2)}});
    
  std::ofstream input = get_input_handle();
  write_rem_section(input, keys);
  write_molecule_section(input);
  input.close();
  exec_qchem();

  /*
    If this is the first call, don't actually read from $QC; there
    will be no overlap to read. It is *required*, however, that the
    above call be made as it saves the relevant overlap quantities for
    the next call.
  */
  if (first_wfoverlap){
    U.eye();
  }
  else{
    readQFMan(QCFILE::FILE_WF_OVERLAP, U);
  }

  first_wfoverlap = false;
}


void QM_QChem::get_ground_energy(arma::vec & e){
  // Build job if we need to
  if (! called(Q::scfman)){
    std::ofstream input = get_input_handle();
    write_rem_section(input, {{"jobtype","sp"}});
    write_molecule_section(input);
    input.close();
    exec_qchem();
  }

  parse_energies(e);
}


// ground + excited states
void QM_QChem::get_all_energies(arma::vec & e){
  if (! called(Q::setman)){
    REMKeys k = excited_rem(false);
    k.insert({
        {"jobtype","sp"},
        {"cis_rlx_dns", "1"}}
      );
    
    std::ofstream input = get_input_handle();
    write_rem_section(input, k);
    write_molecule_section(input);
    input.close();
    called(Q::scfman); // because it will be after this execution
    exec_qchem();
  }

  parse_energies(e);
}


void QM_QChem::get_gradient(arma::mat &g_qm, arma::mat &g_mm, arma::uword surface){
  get_gradient(g_qm, surface);
  if (NMM > 0){
    parse_mm_gradient(g_mm);
  }
}


void QM_QChem::get_gradient(arma::mat &g_qm, arma::uword surface){
  if (surface > excited_states){
    throw std::invalid_argument("Requested surface not computed!");
  }

  std::ofstream input = get_input_handle();
  REMKeys keys = {{"jobtype","force"}};

  if (surface != 0){
    REMKeys ex = excited_rem();
    keys.insert(ex.begin(), ex.end());
    keys.insert({{"cis_state_deriv", std::to_string(surface)}});
  }
  else{
    called(Q::scfman);
  }
  
  write_rem_section(input, keys);  
  write_molecule_section(input);
  input.close();

  exec_qchem();

  readQFMan(QCFILE::FILE_NUCLEAR_GRADIENT, g_qm);
}


/*
  FIXME: Want to pick separate spin states (rem_boys_cis_spin_separate)
  diabatization_method \in {boys, er}
*/
REMKeys QM_QChem::diabatization_rem(std::ofstream & input, size_t I, size_t J){
  //FIXME: need to be careful about NAC skips and excited rem...
  REMKeys keys = excited_rem(false);
  keys.insert({{diabatization_method + "_cis_numstate", "2"},
	       {"calc_nac", "true"},
               {"cis_der_numstate", "2"},  // number of states for DC calculation
	       {"loc_dia_grad", "1"},      // calculates Hamiltonian gradients of diabatic states
	       {"loc_dia_der_type", strictly_diabatic_approximation ? "1": "0"},   // {1 = Strictly diabatic approximation, 0 = exact diabatic dc}
               {"cis_convergence", "7"},
               {"loc_cis_ov_separate", std::to_string(loc_cis_ov_separate)}
    });

  // Necessary for Boys code
  input << "$localized_diabatization" << std::endl;
  input << "(comment)" << std::endl;
  input << I << " " << J << std::endl;
  input << "$end" << std::endl;
  
  input << "$derivative_coupling" << std::endl;
  input << "(comment)" << std::endl;
  input << I << " " << J << std::endl;
  input << "$end" << std::endl;

  return keys;
}


//FIXME: refactor diabatization schemes into their own type
void QM_QChem::get_diabats(arma::cube & gd_qm, arma::mat & U, arma::mat & H, size_t A, size_t B){
  if (spin_diabats.size() > 0){
    get_diabats_spin(gd_qm, U, H);
  }
  else{
    get_diabats_loc(gd_qm, U, H, A, B);
  }
}


void QM_QChem::get_diabats_loc(arma::cube & gd_qm, arma::mat & U, arma::mat & H, size_t I, size_t J){
  if (NMM > 0){
    throw std::logic_error("Diabatic gradients in QM/MM environment not implemented!");
  }
  //if (!called(Q::diabatization)){  // figure out how to make it faster later!
    REMKeys keys = {{"jobtype","sp"}};
    std::ofstream input = get_input_handle();
    REMKeys diab = diabatization_rem(input, I, J);
    keys.insert(diab.begin(), diab.end());
  
    write_rem_section(input, keys);
    write_molecule_section(input);
    input.close();

    exec_qchem();
    //}
  parse_track_diabats(gd_qm, U);

  arma::vec adiabats(1 + excited_states);
  parse_energies(adiabats);
  adiabats -= adiabats(0); // switch to excitations
  
  H = (U.t() * arma::diagmat(adiabats(arma::uvec({I,J}))) * U);
}



void QM_QChem::get_diabats_spin(arma::cube & gd_qm, arma::mat & U, arma::mat & H){
  const int mult_orig = qm_multiplicity;
  const int N = spin_diabats.size();
  U.zeros(N,N);
  H.zeros(N,N);

  bool use_excited_states = diabat_states.size() > 0;
  if (use_excited_states){enable_qink_skips = false;}
  for (int i = 0; i < N; i++){
    qm_multiplicity = spin_diabats[i];
    int state = 0;
    if (use_excited_states){
      state = diabat_states[i];
    }
    
    // when we do the above, will need to make sure we interact with called() correctly

    get_gradient(gd_qm.slice(i), state);

    arma::vec e(excited_states +1, arma::fill::zeros);
    parse_energies(e);

    // Assume zero coupling
    H(i,i) = e(state);
  }

  
  qm_multiplicity = mult_orig;
}


void QM_QChem::parse_track_diabats(arma::cube & gd_qm, arma::mat & U){
  // Read-in diabatic gradients and determine if they need to be swapped
  {
    const size_t NRoot = excited_states;
    const size_t NRoot2 = NRoot * NRoot;
    const size_t NCoord = 3 * NQM;

    // Verified these offsets are correct for the me-bridged closs system
    // c.f. code in $QC/drvman/do_cis_dia_couple.C
    for (int I = 0; I < 2; I++){
      int offset = NCoord*((I+1)*(NRoot + 1) + 4*NRoot2);
      if (strictly_diabatic_approximation){
        offset += 2*NRoot2*NCoord;
      }
      readQFMan(QCFILE::FILE_DERCOUP, gd_qm.slice(I), offset);
      //readQFMan(QCFILE::FILE_DERCOUP, gd_qm.slice(0), (1 +   NRoot + 6*NRoot2)*NCoord);
      //readQFMan(QCFILE::FILE_DERCOUP, gd_qm.slice(1), (2 + 2*NRoot + 6*NRoot2)*NCoord);
    }
    
    readQFMan(QCFILE::FILE_DIAB_ROT_MAT, U);
  }

  /*
    We track the identity of the diabatic gradients (states) by
    tracking the dipole centers associated with each state. On each
    invocation, we construct all inter-state distances, (D)_{A',B} =
    ||{\mu}_{A'A'} - {\mu}_{BB}|| (where the primed indicies are the
    previous (reference) dipole centers) and taking as identical the
    state with the smallest difference. The other state is therefore
    also determined.

    This information is not exposed, but all gradients returned by
    repeated invocation of this function are consistent up to the
    criterion above.
  */
  
  size_t numkeep2 = 2*2; // keeping 2 diabatic states
  arma::mat mu(numkeep2,3);
  // our little corner of FILE_DC_DIPS is laid out as:
  // u*(N*N) + A*N + B; u \on {x,y,z}; A,B \on 1..N
  // numkeep2 = N*N = 4
  // yields (Mu)_{ABu}
  readQFMan(QCFILE::FILE_DC_DIPS, mu, 3 + 9 * numkeep2);
  mu = mu.t();  // now each column is a dipole: AA, AB, BA, BB
  
  mu.save(arma::hdf5_name("gifs.hdf5",
                          "/diabaticdipoles/" + std::to_string(call_idx()),
                          arma::hdf5_opts::replace));
  
  mu = mu.cols(arma::uvec({0,3}));  // now the columns are AA BB

  if (!(donor_acceptor_ref.n_cols == 2 && donor_acceptor_ref.n_rows == 3)){
    std::cerr << "[QM_QChem] " << call_idx() << ": No reference dipoles provided via 'donor_acceptor_ref'; this trajectory will not be reproducible." << std::endl;
    donor_acceptor_ref = mu;
  }

  arma::mat D(2,2);
  
  for (arma::uword i = 0; i < 2; i++){
    for (arma::uword j = 0; j < 2; j++){
      D(i,j) = arma::norm(donor_acceptor_ref.col(i)-mu.col(j));
    }
  }

  size_t idx=arma::index_min(arma::vectorise(D));
  if ((1==idx) || (2==idx)){ // if the smallest difference is off-diagonal
    gd_qm.slice(0).swap(gd_qm.slice(1));
    U.swap_rows(0,1);
    mu.swap_cols(0,1);
    std::cerr << "[QM_QChem] " << call_idx() << ": Swapping diabats; " << D(idx) << " < " << D(0) << std::endl;
  }

  donor_acceptor_ref = mu; // update for subsequent invocations
}


void QM_QChem::get_nac_vector(arma::mat & nac, size_t I, size_t J){
  /*
    cannot skip setman for nac; no need to skip scfman because it will
    be so quick
  */
  REMKeys keys = excited_rem(false);

  /*
    relaxed density should be computed automatically, but we set it
    explicitly here since we're setting scf- and setman as called()
  */
  called(Q::scfman); called(Q::setman);
  keys.insert({
      {"jobtype","sp"},
      {"calc_nac", "true"},
      {"cis_der_numstate", "2"}, // Always, in our case, between 2 states
      {"cis_rlx_dns", "1"}
    });

  if (NMM > 0){
    keys.insert({{"nac_pointcharge", "1"},
                 {"qm_mm", "true"}});
  }
  
  std::ofstream input = get_input_handle();
  write_rem_section(input, keys);

  /*
    Section format:
    $derivative_coupling
      comment line
      A B ...
    $end
    where A, B, ... are the states between which to compute the
    derivative coupling. The $rem section must include
    cis_der_numstates equal to the total number of states
    specified. The ground state is 0.

  */
  input << "$derivative_coupling" << std::endl;
  input << "comment" << std::endl;
  input << I << " " << J << std::endl;
  input << "$end" << std::endl;
  
  write_molecule_section(input);
  input.close();

  exec_qchem();

  readQFMan(QCFILE::FILE_DERCOUP, nac);
}


void QM_QChem::parse_energies(arma::vec &e){
  readQFMan(QCFILE::FILE_ENERGY, e.memptr(), 1, POS_CRNT_TOTAL_ENERGY);

  if (excited_states > 0 && e.n_elem > 1){
    double e_ground = e[0];

    if (!(e.n_elem >= excited_states + 1)){
      throw std::range_error("Insufficient space for all excited states in vec e!");
    }

    size_t count = readQFMan(QCFILE::FILE_SET_ENERGY, e.memptr() + 1, excited_states, POS_BEGIN);
    if (count != excited_states){
      throw std::runtime_error("Unable to parse energies!");
    }

    // Energies for higher states given in terms of excitations so we
    // must add in the ground.
    e += e_ground;
    e[0] = e_ground;
  }
  e.save(arma::hdf5_name("gifs.hdf5",
                         "/energies/" + std::to_string(call_idx()),
                         arma::hdf5_opts::replace));
}


void QM_QChem::parse_mm_gradient(arma::mat &g_mm){
  /* Size of g_mm updated above in GifsImpl::update_gradient() */
  readQFMan(QCFILE::FILE_EFIELD, g_mm);

  for (size_t i = 0; i < NMM; i++){
    const double q = chg_mm[i];
    // FILE_EFIELD really contains the *field* so we need a negative to get the gradient.
    g_mm.col(i) *= -1.0 * q;
  }
  /*
    FIXME: prevent efield.dat from being written with REM_QMMM_EXT_GIFS set

    This method replaced parsing efield.dat, which has the following
    format:

    Ex(mm1) Ey(mm1) Ez(mm1)
    ...
    Ex(mmN) Ey(mmN) Ez(mmN)
    Ex(qm1) Ey(qm1) Ez(qm1)
    ...
    Ex(qmN) Ey(qmN) Ez(qmN)

    where Ea(u) is the component of the electric field in the a
    direction at u and mmi and qmi are the coordinates of the ith mm
    and qm atoms respectively. N.N. The product of the field and the
    charge is a force rather than the gradient.
  */
}


const std::string QM_QChem::get_qcprog(void){
  char * qc_str = std::getenv("QC");
  std::string default_path = "qcprog.exe";
  if (nullptr == qc_str){
    std::cerr << "Warning, $QC not set; using " << default_path << std::endl;
    return default_path;
  }
  else{
    return std::string(qc_str) + "/exe/qcprog.exe";
  }
}


const std::string QM_QChem::get_qcwdir(void){
  char * pwd = std::getenv("PWD");
  std::string qcwdir = std::string(pwd) + "/GQSH";

  static bool create = true;

  if (create){
    struct stat st = {};
    if (-1 == stat(qcwdir.c_str(), &st)){
      mkdir(qcwdir.c_str(), 0700);
    }
    create = false;
  }

  return qcwdir;
}


/*
  In order of preference:
  1)               If $QCSCRATCH is set and an absolute path use it
  2) Failing that, if conf_dir   is set and an absolute path use it
  3) Failing that, use our default path

  The above are set in reverse.
*/
const std::string QM_QChem::get_qcscratch(std::string conf_dir){
  char * pwd = std::getenv("PWD");
  std::string scratch_path = std::string(pwd) + "/GQSH.sav";

  // if the config directory is valid, take it as the default
  if (conf_dir[0] == '/'){
    scratch_path = conf_dir;
  }
  
  char * qc_str = std::getenv("QCSCRATCH");
  if (nullptr == qc_str){
    std::cerr << "Warning, $QCSCRATCH not set; using " << scratch_path << std::endl;
  }
  else if (qc_str[0] != '/'){
    std::cerr << "Warning, $QCSCRATCH is not an absolute path; instead using " << scratch_path << std::endl;
  }
  else{
    scratch_path = std::string(qc_str);
    // If using a global scratch path, make it unique
    // FIXME: could still collide if QCSCRATCH is on a network device
    std::srand(getpid());
    scratch_path += "/INAQS_" + std::to_string(std::time(nullptr)) + "_" + std::to_string(rand());
  }

  /* Make sure the directory exists; create otherwise */
  struct stat st = {};
  if (-1 == stat(scratch_path.c_str(), &st)){
    if (mkdir(scratch_path.c_str(), 0700)){
      throw std::runtime_error("Unable to construct scratch directory at " + scratch_path);
    }
  }

  /*
    Set $QCSCRATCH to whatever we're using. Incidentally, this fixes a
    corner case where Q-Chem will crash if QCSCRATCH=""; set, but to
    the empty string.
  */
  setenv("QCSCRATCH", scratch_path.c_str(), true);
  std::cerr << "taking scratch path to be: " << scratch_path << std::endl;

  return scratch_path;
}

void QM_QChem::exec_qchem(void){
  if (dump_qc_output){
    qc_log_file = qc_log_file_idx();
  }

  std::string cmd = "cd " + get_qcwdir() + "; " + // change to  target WD
    qc_executable + " " +
    qc_scratch_directory + "/" + qc_input_file + " " +
    qc_scratch_directory + " >" + qc_log_file + " 2>&1";
  int status = std::system(cmd.c_str());
  if (status){
    throw std::runtime_error("Q-Chem could not be called or exited abnormally; see " + get_qcwdir() + "/" + qc_log_file);
  }
  first_call = false;
}


/* Consistent formatting and file name */
std::ofstream QM_QChem::get_input_handle(){
  std::ofstream os(qc_scratch_directory + "/" + qc_input_file);
  os.setf(std::ios_base::fixed, std::ios_base::floatfield);
  os.precision(std::numeric_limits<double>::digits10);

  if (!os){
    throw std::runtime_error("Unable to open for writing the file: " + qc_scratch_directory + "/" + qc_input_file);
  }
  
  return os;
}

/*
  detectQinks must be false if in constructing your job and before
  calling exec_qchem() you:
  - call called() at all
  - call excited_rem() more than once
*/
REMKeys QM_QChem::excited_rem(bool detectQinks){
  REMKeys keys
    {
      {"cis_n_roots", std::to_string(excited_states)},
      {"set_roots_orig", std::to_string(excited_states)}, // make sure we always use the same number of states
      {"max_cis_cycles", "500"}, // same as set_iter
      {"cis_singlets", std::to_string(singlets)},  
      {"cis_triplets", std::to_string(triplets)} 
    };

  if (detectQinks && enable_qink_skips){
    bool skip_scfman = called(Q::scfman);  // called() updates internal state so we
    bool skip_setman = called(Q::setman);  // need to make sure both are touched.
    if (skip_scfman){
      keys.insert({{"skip_scfman", "1"}});
      if (skip_setman){
        keys.insert({{"skip_setman", "1"}});
      }
    }
    else{
      keys.insert({{"cis_rlx_dns", "1"}});
    }
  }

  return keys;
}


void QM_QChem::write_rem_section(std::ostream &os, const REMKeys &options){
  // Default options 
  REMKeys rem_keys
    {
      {"method",                    exchange_method},
      {"basis",                     basis_set},
      {"scf_algorithm",             scf_algorithm},
      {"scf_convergence",           std::to_string(scf_convergence)},
      {"thresh",                    std::to_string(two_electron_thresh)},
      {"thresh_diis_switch",        "7"},     // control change-over to GDM if
      {"max_diis_cycles",           "50"},    // scf_alg is DIIS_GDM
      {"sym_ignore",                "true"},
      {"max_scf_cycles",            "500"},
      {"unrestricted",              std::to_string(unrestricted)},
      {"cis_s2_thresh",             std::to_string(int(sing_thresh * 100))},
      {"print_input",               std::to_string(dump_qc_input)},
      {"input_bohr",                "true"}  // works for MM charges; c.f. $QC/libgen/PointCharges.C
    };

  if (NMM > 0){
    rem_keys.emplace("qm_mm", "true");
    rem_keys.emplace("skip_charge_self_interact", "1");  // since the MD driver will compute MM-MM interaction
  }
    
  if (spin_flip){
    rem_keys.emplace("spin_flip", "1");
    rem_keys.emplace("sts_mom", "1");
  }
  if (first_call){
    rem_keys.emplace("qmmm_ext_gifs", "1");
  }
  else{
    rem_keys.emplace("qmmm_ext_gifs", "2");
    rem_keys.emplace("scf_guess",     scf_guess);
  }  

  /*
    Make sure passed options take precedence?  Could do this by
    inserting the defaults into the passed options rather than the
    other way around.
  */
  
  // merge in options; usually will include jobtype 
  rem_keys.insert(options.begin(), options.end());
  
  os << "$rem" << std::endl;
  for (const auto& e: rem_keys){
    os << e.first << " " << e.second << std::endl;
  }
  os << "$end" << std::endl;

  /*
    write the contents of the verbatim file.

    N.B. that if the verbatim file contains a $rem section, those
    options will take precedence over those written above.
  */
  if (verbatim_file.length() > 0)
  {
    std::ifstream verb(verbatim_file);
    if (!verb){
      throw std::runtime_error("Unable to open verbatim file for reading: " + verbatim_file);
    }
    std::string buff;
    while (std::getline(verb, buff)) {
      os << buff << std::endl;
    }
    verb.close();
  }
}


/*
  FIXME: should do better formatting with std::right << std::fixed <<
  std::setprecision(precision) << std::setw(precision+5) and friends
*/
void QM_QChem::write_molecule_section(std::ostream &os){
  /*
    format of $molecule & $externall_charges sections:
    $molecule
    [charge] [multiplicity]
    [atomic-number] [x-coord] [y-coord] [x-coord]
    ...
    $end

    $external_charges
    [x-coord] [y-coord] [x-coord] [charge]
    ...
    $end
  */  
  os << "$molecule" << std::endl;
  os << qm_charge << " " << qm_multiplicity << std::endl;

  for (size_t i = 0; i < NQM; i++){
    os << atomids[i]      << " ";        // id
    os << crd_qm[i*3 + 0] << " ";        // x
    os << crd_qm[i*3 + 1] << " ";        // y
    os << crd_qm[i*3 + 2] << std::endl;  // x
  }
  os <<  "$end" << std::endl;

  if (NMM > 0){
    if(externalcharges_hack){
      std::string extfname = "external_charges.in";
      os << std::endl << "$external_charges" << std::endl;
      os << extfname << std::endl;
      os << "$end" << std::endl;

      std::ofstream osext(get_qcwdir() + "/" + extfname);
      osext.setf(std::ios_base::fixed, std::ios_base::floatfield);
      osext.precision(std::numeric_limits<double>::digits10);
      write_external_charges_section(osext);
      osext.close();
    }
    else{
      write_external_charges_section(os);
    }
  }
}

void QM_QChem::write_external_charges_section(std::ostream &os){
  os << std::endl << "$external_charges" << std::endl;
  for (size_t i = 0; i < NMM; i++){
    os << crd_mm[i*3 + 0] << " ";        // x
    os << crd_mm[i*3 + 1] << " ";        // y
    os << crd_mm[i*3 + 2] << " ";        // z
    os << chg_mm[i]       << std::endl;  // charge
  }
  os << "$end" << std::endl;
}

/*
  Sentinel system to track whether the respective q-chem qink has been
  called since the last call to update(). See/update the protected
  nested enum class Q (as in sentinel) for the list of relevant qinks
*/
bool QM_QChem::called(Q q){
  static std::unordered_map<Q, int> calls;
  if(call_idx() == calls[q]){
    return true;
  }
  else{
    calls[q] = call_idx();
    return false;
  }
}


/*
  Given a q-qchem file number (see qm_qchem.hpp for examples), read N
  elements from the file (in the scratch directory) starting from the
  offset into memptr, which must point to a block of memory with space
  for at least N doubles.
*/
size_t QM_QChem::readQFMan(QCFILE filenum, double * memptr, size_t N, size_t offset){
  std::string path = qc_scratch_directory + "/" + std::to_string((int) filenum) + ".0";
  std::ifstream ifile;
  ifile.open(path, std::ios::in | std::ios::binary);

  if (!ifile.is_open()){
    throw std::ios_base::failure("Error: cannot read from " + path);
  }

  char buffer[sizeof(double)];
  size_t i = 0;
  size_t count = 0;
  while(ifile.read(buffer, sizeof(double))){
    if ((i >= offset) && (i < offset + N)){
      double * d = (double *) buffer;
      *(memptr + count++) = *d;
    }
    i++;
  }
  ifile.close();

  return count;
}

std::string QM_QChem::to_string(const QCFILE f){
  std::string name;
  switch (f){
  case QCFILE::FILE_SET_ENERGY:
    name = "FILE_SET_ENERGY, 72.0";
    break;
  case QCFILE::FILE_ENERGY:
    name = "FILE_ENERGY, 99.0";
    break;
  case QCFILE::FILE_NUCLEAR_GRADIENT:
    name = "FILE_NUCLEAR_GRADIENT, 131.0";
    break;
  case QCFILE::FILE_EFIELD:
    name = "FILE_EFIELD, 329.0";
    break;
  case QCFILE::FILE_WF_OVERLAP:
    name = "FILE_WF_OVERLAP, 398.0";
    break;
  case QCFILE::FILE_DIAB_ROT_MAT:
    name = "FILE_DIAB_ROT_MAT, 941.0";
    break;
  case QCFILE::FILE_TRANS_DIP_MOM:
    name = "FILE_TRANS_DIP_MOM, 942.0";
    break;
  case QCFILE::FILE_DC_DIPS:
    name = "FILE_DC_DIPS, 966.0";
    break;
  case QCFILE::FILE_DERCOUP:
    name = "FILE_DERCOUP, 967.0";
    break;
  case QCFILE::FILE_CIS_S2:
    name = "FILE_CIS_S2, 1200.0";
    break;
  }
  return name;
}
