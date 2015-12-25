#include <kvs/Message>
#include <kvs/PolygonObject>
#include <kvs/glut/Application>
#include <kvs/glut/Screen>
#include <kvs/UnstructuredVolumeObject>
#include <kvs/Endian>
#include <kvs/PointObject>
#include <kvs/CellByCellMetropolisSampling>
#include <kvs/ParticleBasedRenderer>
#include <kvs/TransferFunction>
#include <kvs/Vector3>
#include <kvs/StochasticPolygonRenderer>
#include <kvs/StochasticTetrahedraRenderer>
#include <kvs/StochasticRenderingCompositor>
#include <kvs/Scene>
#include <kvs/ObjectManager>
#include <kvs/RendererManager>
#include <iostream>
#include <fstream>

//半径方向の分割数(ocore)
const size_t Rad_n = 201;

//緯度方向の分割数(ocore)
const size_t Lat_n = 204;

//経度方向の分割数(ocore)
const size_t Lon_n = 608;

//一辺の分割数(icore)
const size_t Zhong_n = 222;

//数値データのベクトル長
const size_t veclen = 1;

//接点の数(ocore)
const size_t nnodes_o = Rad_n * Lat_n * Lon_n;

//要素数(ocore)
const size_t ncells_o = (Rad_n - 1) * (Lat_n - 1) * (Lon_n - 1);

//接点の数(icore)
const size_t nnodes_i = Zhong_n * Zhong_n * Zhong_n;

//要素数(icore)
const size_t ncells_i = (Zhong_n - 1) * (Zhong_n - 1) * (Zhong_n - 1);

//円周率
const float Pi = 3.141593f;

//r(ocore)
const float R_max = 1.0f;
const float R_min = 0.35f;
const float R_d   = ( R_max - R_min ) / ( Rad_n - 1 );

//θ(ocore)
const float Th_max = Pi - Pi / 4.0f;
const float Th_min = Pi / 4.0f;
const float Th_d   = ( Th_max - Th_min ) / ( Lat_n -3 );

//φ(ocore)
const float Ph_max = ( 3 * Pi ) / 4.0f;
const float Ph_min = -( 3 * Pi ) / 4.0f;
const float Ph_d   = ( Ph_max - Ph_min ) / ( Lon_n - 5 );

//ix(= iy, iz),dix(= diy, diz)
const float R_i = R_min + R_d * 2;
const float ix_max = R_i;
const float ix_min = -R_i;
const float dix = ( ix_max - ix_min ) / ( Zhong_n - 1 );

//---SetYinCoordArray-----------------------------------------------------------------
kvs::Real32* SetYinCoordArray()
{
  //座標値配列を格納する配列
  static kvs::Real32* CoordArray = new kvs::Real32 [nnodes_o * 3];

  //coordArrayの配列の要素番号
  static int coord_i = 0;

  //グリッド交差点の座標決定
  for(int k = 0; k < (int)Lon_n; k++){
    float ph = Ph_min + Ph_d * (float)(k - 2); 
    for(int j = 0; j < (int)Lat_n; j++){
      float th = Th_min + Th_d * (float)(j - 1);
      for(int i = 0; i < (int)Rad_n; i++){
	float r = R_min + R_d * (float)i;
	float x = r * sinf(th) * cosf(ph);
	float y = r * sinf(th) * sinf(ph);
	float z = r * cosf(th);
	CoordArray[coord_i] = x;
	CoordArray[coord_i + 1] = y;
	CoordArray[coord_i + 2] = z;
	coord_i += 3;
      }
    }
  }
  //座標値配列を返す
  return(CoordArray);
}  
//---SetYinCoordArray----------------------------------------------------------------

//---SetYangCoordArray---------------------------------------------------------------
kvs::Real32* SetYangCoordArray()
{
  //座標値配列を格納する配列
  static kvs::Real32* CoordArray = new kvs::Real32 [nnodes_o * 3];

  //coordArrayの配列の要素番号
  static int coord_i = 0;

  //グリッド交差点の座標決定
  for(int k = 0; k < (int)Lon_n; k++){
    float ph = Ph_min + Ph_d * (float)(k - 2);
    for(int j = 0; j < (int)Lat_n; j++){
      float th = Th_min + Th_d * (float)(j - 1);
      for(int i = 0; i < (int)Rad_n; i++){
        float r = R_min + R_d * (float)i;
	float x = r * sinf(th) * cosf(ph);
	float y = r * sinf(th) * sinf(ph);
	float z = r * cosf(th);
        CoordArray[coord_i] = -x;
        CoordArray[coord_i + 1] = z;
        CoordArray[coord_i + 2] = y;
        coord_i += 3;
      }
    }
  }
  //座標値配列を返す                               
  return(CoordArray);
}
//---SetYangCoordArray------------------------------------------------------------

//---SetZhongCoordArray-----------------------------------------------------------------
kvs::Real32* SetZhongCoordArray()
{
  //座標値配列を格納する配列
  kvs::Real32* CoordArray = new kvs::Real32 [nnodes_i * 3];

  //coordArrayの配列の要素番号
  static int coord_i = 0;

  //グリッド交差点の座標決定
  for(int k = 0; k < (int)Zhong_n; k++){
    for(int j = 0; j < (int)Zhong_n; j++){
      for(int i = 0; i < (int)Zhong_n; i++){
	float x = ix_min + dix * i;
	float y = ix_min + dix * j;
	float z = ix_min + dix * k;
	CoordArray[coord_i] = x;
	CoordArray[coord_i + 1] = y;
	CoordArray[coord_i + 2] = z;
	coord_i += 3;
      }
    }
  }
  //座標値配列を返す
  return(CoordArray);
}  
//---SetZhongCoordArray----------------------------------------------------------------

//---SetYinValueArray----------------------------------------------------------------
//バイナリファイルからvalueに格納
float* SetYinValueArray(){
  // 数値データを格納する配列
  float* ValueArray = new float [nnodes_o * veclen];

  //バイナリファイルの読み取り
  std::ifstream ifs("../bx_vx/oct09b.011.wyin.vx.n000250000.t00302", std::ios_base::in | std::ios_base::binary);

  //ファイル読み取りエラー処理
  if(!ifs){
    std::cerr << "バイナリファイルを開けませんでした。" << std::endl;
    exit(-1);
  }

  //数値データの格納
  ifs.seekg(4,std::ios::beg);
  ifs.read((char*)ValueArray, nnodes_o * veclen * sizeof(float));

  //エンディアン
  kvs::Endian::Swap(ValueArray, nnodes_o * veclen);

  //数値データの配列を返す
  return(ValueArray);
}
//---SetYinValueArray---------------------------------------------------------------

//---SetYangValueArray----------------------------------------------------------------
//バイナリファイルからvalueに格納
float* SetYangValueArray(){

  // 数値データを格納する配列
  float* ValueArray = new float [nnodes_o * veclen];

  //バイナリファイルの読み取り
  std::ifstream ifs("../bx_vx/oct09b.011.wyng.vx.n000250000.t00302", std::ios_base::in | std::ios_base::binary);

  //ファイル読み取りエラー処理
  if(!ifs){
    std::cerr << "バイナリファイルを開けませんでした。" << std::endl;
    exit(-1);
  }

  //数値データの格納
  ifs.seekg(4,std::ios::beg);
  ifs.read((char*)ValueArray, nnodes_o * veclen * sizeof(float));

  //エンディアン
  kvs::Endian::Swap(ValueArray, nnodes_o * veclen);

  //数値データの配列を返す
  return(ValueArray);
}
//---SetYangValueArray---------------------------------------------------------------

//---SetZhongValueArray----------------------------------------------------------------
//バイナリファイルからvalueに格納
float* SetZhongValueArray(){
  // 数値データを格納する配列
  float* ValueArray = new float [nnodes_i * veclen];

  //バイナリファイルの読み取り
  std::ifstream ifs("../bx_vx/oct09b.011.icore_3d.vx.n000250000.t00302", std::ios_base::in | std::ios_base::binary);

  //ファイル読み取りエラー処理
  if(!ifs){
    std::cerr << "バイナリファイルを開けませんでした。" << std::endl;
    exit(-1);
  }

  //数値データの格納
  ifs.seekg(4,std::ios::beg);
  ifs.read((char*)ValueArray, nnodes_i * veclen * sizeof(float));

  //エンディアン
  kvs::Endian::Swap(ValueArray, nnodes_i * veclen);

  //数値データの配列を返す
  return(ValueArray);
}
//---SetZhongValueArray---------------------------------------------------------------

//---SetYinYangConnectionArray--------------------------------------------------------------
kvs::UInt32* SetYinYangConnectionArray(){
  //接続情報を格納する配列
  kvs::UInt32* ConnectionArray = new kvs::UInt32 [ncells_o * 8];
  //ConnectionArrayの要素番号
  int connect_i = 0;

  //指定ノードの番号
  int node_n = 0;

  //接続情報
  for(int k = 0; k < (int)Lon_n - 1; k++){
    for(int j = 0; j < (int)Lat_n - 1; j++){
      for(int i = 0; i < (int)Rad_n - 1; i++){
	ConnectionArray[connect_i] = node_n;
	ConnectionArray[connect_i + 1] = node_n + 1;
	ConnectionArray[connect_i + 2] = node_n + (Rad_n * Lat_n) + 1;
	ConnectionArray[connect_i + 3] = node_n + (Rad_n * Lat_n);
	ConnectionArray[connect_i + 4] = node_n + Rad_n;
	ConnectionArray[connect_i + 5] = node_n + Rad_n + 1;
	ConnectionArray[connect_i + 6] = node_n + Rad_n + (Rad_n * Lat_n) + 1;
	ConnectionArray[connect_i + 7] = node_n + Rad_n + (Rad_n * Lat_n);
	connect_i += 8;
	node_n++;
      }
      node_n++;
    }
    node_n += (int)Rad_n;
  }

  return(ConnectionArray);
}
//---SetYinYangConnectionArray----------------------------------------------------------

//---SetZhongConnectionArray----------------------------------------------------------
kvs::UInt32* SetZhongConnectionArray(){
  //接続情報を格納する配列
  kvs::UInt32* ConnectionArray = new kvs::UInt32 [ncells_i * 8];

  //ConnectionArrayの要素番号
  static int connect_i = 0;

  //指定ノードの番号
  static int node_n = 0;

  //接続情報
  for(int k = 0; k < (int)Zhong_n - 1; k++){
    for(int j = 0; j < (int)Zhong_n - 1; j++){
      for(int i = 0; i < (int)Zhong_n - 1; i++){
	ConnectionArray[connect_i] = node_n;
	ConnectionArray[connect_i + 1] = node_n + 1;
	ConnectionArray[connect_i + 2] = node_n + (Zhong_n * Zhong_n) + 1;
	ConnectionArray[connect_i + 3] = node_n + (Zhong_n * Zhong_n);
	ConnectionArray[connect_i + 4] = node_n + Zhong_n;
	ConnectionArray[connect_i + 5] = node_n + Zhong_n + 1;
	ConnectionArray[connect_i + 6] = node_n + Zhong_n + (Zhong_n * Zhong_n) + 1;
	ConnectionArray[connect_i + 7] = node_n + Zhong_n + (Zhong_n * Zhong_n);
	connect_i += 8;
	node_n++;
      }
      node_n++;
    }
    node_n += (int)Zhong_n;
  }

  return(ConnectionArray);
}
//---SetZhongConnectionArray----------------------------------------------------------

//---CreateYinUnstructureVolumeObject-----------------------------------------------
// 非構造格子型ボリュームオブジェクト（四面体１次要素）を生成する関数
kvs::UnstructuredVolumeObject* CreateYinUnstructuredVolumeObject()
{
  // KVSの配列クラス（kvs::ValueArray）および動的型配列クラス（kvs::AnyValueArray）にセットする。
  kvs::ValueArray<kvs::Real32> coords( SetYinCoordArray(), nnodes_o * 3);
  kvs::ValueArray<kvs::UInt32> connections( SetYinYangConnectionArray(), ncells_o * 8);
  kvs::ValueArray<float> values( SetYinValueArray(), nnodes_o * veclen);

  // 非構造格子型ボリュームオブジェクトを生成する。
  kvs::UnstructuredVolumeObject* object = new kvs::UnstructuredVolumeObject();
  object->setCellType( kvs::UnstructuredVolumeObject::Hexahedra );
  object->setVeclen( veclen );
  object->setNumberOfNodes( nnodes_o );
  object->setNumberOfCells( ncells_o );
  object->setCoords( coords );
  object->setConnections( connections );
  object->setValues( values );
  object->updateMinMaxValues();

  return object;
}
//---CreateYinUnstructureVolumeObject-----------------------------------------------

//---CreateYangUnstructureVolumeObject-----------------------------------------------
// 非構造格子型ボリュームオブジェクト（四面体１次要素）を生成する関数
kvs::UnstructuredVolumeObject* CreateYangUnstructuredVolumeObject()
{
  // KVSの配列クラス（kvs::ValueArray）および動的型配列クラス（kvs::AnyValueArray）にセットする。
  kvs::ValueArray<kvs::Real32> coords( SetYangCoordArray(), nnodes_o * 3);
  kvs::ValueArray<kvs::UInt32> connections( SetYinYangConnectionArray(), ncells_o * 8);
  kvs::ValueArray<float> values( SetYangValueArray(), nnodes_o * veclen);

  // 非構造格子型ボリュームオブジェクトを生成する。
  kvs::UnstructuredVolumeObject* object = new kvs::UnstructuredVolumeObject();
  object->setCellType( kvs::UnstructuredVolumeObject::Hexahedra );
  object->setVeclen( veclen );
  object->setNumberOfNodes( nnodes_o );
  object->setNumberOfCells( ncells_o );
  object->setCoords( coords );
  object->setConnections( connections );
  object->setValues( values );
  object->updateMinMaxValues();

  return object;
}
//---CreateYangUnstructureVolumeObject-----------------------------------------------

//---CreateZhongUnstructureVolumeObject-----------------------------------------------
// 非構造格子型ボリュームオブジェクト（四面体１次要素）を生成する関数
kvs::UnstructuredVolumeObject* CreateZhongUnstructuredVolumeObject()
{
  // KVSの配列クラス（kvs::ValueArray）および動的型配列クラス（kvs::AnyValueArray）にセットする。
  kvs::ValueArray<kvs::Real32> coords( SetZhongCoordArray(), nnodes_i * 3);
  kvs::ValueArray<kvs::UInt32> connections( SetZhongConnectionArray(), ncells_i * 8);
  kvs::ValueArray<float> values( SetZhongValueArray(), nnodes_i * veclen);

  // 非構造格子型ボリュームオブジェクトを生成する。
  kvs::UnstructuredVolumeObject* object = new kvs::UnstructuredVolumeObject();
  object->setCellType( kvs::UnstructuredVolumeObject::Hexahedra );
  object->setVeclen( veclen );
  object->setNumberOfNodes( nnodes_i );
  object->setNumberOfCells( ncells_i );
  object->setCoords( coords );
  object->setConnections( connections );
  object->setValues( values );
  object->updateMinMaxValues();

  return object;
}
//---CreateZhongUnstructureVolumeObject-----------------------------------------------

//---main------------------------------------------------------------------------
// メイン関数
int main( int argc, char** argv )
{
  kvs::glut::Application app( argc, argv );

  // 非構造格子型ボリュームオブジェクト（四面体１次要素）の生成
  kvs::UnstructuredVolumeObject* volumeYin = CreateYinUnstructuredVolumeObject();
  kvs::UnstructuredVolumeObject* volumeYang = CreateYangUnstructuredVolumeObject();
  kvs::UnstructuredVolumeObject* volumeZhong = CreateZhongUnstructuredVolumeObject();

  //数値の最大最小の再設定 
  double min_yyz  = kvs::Math::Min( volumeYin->minValue(), volumeYang->minValue(), volumeZhong->minValue() );
  double max_yyz  = kvs::Math::Max( volumeYin->maxValue(), volumeYang->maxValue(), volumeZhong->maxValue() );

  //valueの最大最小の設定
  volumeYin->setMinMaxValues( min_yyz, max_yyz );
  volumeYang->setMinMaxValues( min_yyz, max_yyz );
  volumeZhong->setMinMaxValues( min_yyz, max_yyz );


  //領域の最大最小の設定
  volumeYin->updateMinMaxCoords();
  volumeYang->updateMinMaxCoords();
  volumeZhong->updateMinMaxCoords();

  std::cout << "( min_yin, min_yang, min_zhong ) = ( " << volumeYin->minObjectCoord() << ", " << volumeYang->minObjectCoord() << ", " << volumeZhong->minObjectCoord() << " )" << std::endl;

  const kvs::Vec3 min_yin = volumeYin->minObjectCoord();
  const kvs::Vec3 min_yang = volumeYang->minObjectCoord();
  const kvs::Vec3 min_zhong = volumeZhong->minObjectCoord();

  const kvs::Vec3 max_yin = volumeYin->maxObjectCoord();
  const kvs::Vec3 max_yang = volumeYang->maxObjectCoord();
  const kvs::Vec3 max_zhong = volumeZhong->maxObjectCoord();

  const float min_x = kvs::Math::Min( min_yin.x(), min_yang.x(), min_zhong.x() );
  const float min_y = kvs::Math::Min( min_yin.y(), min_yang.y(), min_zhong.y() );
  const float min_z = kvs::Math::Min( min_yin.z(), min_yang.z(), min_zhong.z() );
  const kvs::Vec3 min = kvs::Vec3( min_x, min_y, min_z );

  const float max_x = kvs::Math::Max( max_yin.x(), max_yang.x(), max_zhong.x() );
  const float max_y = kvs::Math::Max( max_yin.y(), max_yang.y(), max_zhong.y() );
  const float max_z = kvs::Math::Max( max_yin.z(), max_yang.z(), max_zhong.z() );
  const kvs::Vec3 max = kvs::Vec3( max_x, max_y, max_z );

  volumeYin->setMinMaxObjectCoords( min, max );
  volumeYang->setMinMaxObjectCoords( min, max );
  volumeZhong->setMinMaxObjectCoords( min, max );

  volumeYin->setMinMaxExternalCoords( min, max );
  volumeYang->setMinMaxExternalCoords( min, max );
  volumeZhong->setMinMaxExternalCoords( min, max );

  //ParticleBasedRendere
  const size_t repetitions = 50;
  const size_t subpixels = 1; // fixed to '1'
  const size_t level = static_cast<size_t>( subpixels * std::sqrt( double( repetitions ) ) );
  const float step = 0.1f;

  //不透明度の設定
  kvs::OpacityMap omap( 256 );
  omap.addPoint( 0, 1.0 );
  omap.addPoint( 90, 0.0 );
  omap.addPoint( 180, 0.0 );
  omap.addPoint( 255, 1.0 );
  omap.create();

  const kvs::TransferFunction tfunc( omap ); //伝達関数 

  kvs::PointObject* objectYin = new kvs::CellByCellMetropolisSampling( volumeYin, level, step, tfunc );
  objectYin->print( std::cout << std::endl );
  delete volumeYin;

  kvs::PointObject* objectYang = new kvs::CellByCellMetropolisSampling( volumeYang, level, step, tfunc );
  objectYang->print( std::cout << std::endl );
  delete volumeYang;

  kvs::PointObject* objectZhong = new kvs::CellByCellMetropolisSampling( volumeZhong, level, step, tfunc );
  objectZhong->print( std::cout << std::endl );
  delete volumeZhong;

  //レンダラーの設定
  kvs::glsl::ParticleBasedRenderer* rendererYin = new kvs::glsl::ParticleBasedRenderer();
  rendererYin->disableShading();

  kvs::glsl::ParticleBasedRenderer* rendererYang = new kvs::glsl::ParticleBasedRenderer();
  rendererYang->disableShading();

  kvs::glsl::ParticleBasedRenderer* rendererZhong = new kvs::glsl::ParticleBasedRenderer();
  rendererZhong->disableShading();

  kvs::glut::Screen screen( &app );
  screen.setTitle( "ParticleBasedRenderer" );
  screen.registerObject( objectYin, rendererYin );
  screen.registerObject( objectYang, rendererYang );
  screen.registerObject( objectZhong, rendererZhong );
  screen.show();

  kvs::StochasticRenderingCompositor compositor( screen.scene() );
  compositor.setRepetitionLevel( repetitions);
  compositor.enableLODControl();
  screen.setEvent( &compositor );

  return app.run();
}
//---main------------------------------------------------------------------------
