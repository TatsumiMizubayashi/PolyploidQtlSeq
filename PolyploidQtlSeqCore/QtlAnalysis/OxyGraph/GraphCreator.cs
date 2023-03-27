using OxyPlot;
using OxyPlot.SkiaSharp;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフクリエーターの抽象クラス
    /// </summary>
    internal abstract class GraphCreator
    {
        /// <summary>
        /// グラフ設定
        /// </summary>
        protected readonly GraphConfig _config;

        /// <summary>
        /// グラフタイトル
        /// </summary>
        private readonly string _title;

        /// <summary>
        /// Y軸タイトル
        /// </summary>
        private readonly string _yAxisTitle;

        /// <summary>
        /// 画像ファイル名
        /// </summary>
        private readonly string _fileBaseName;

        /// <summary>
        /// グラフクリエーターインスタンスを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        /// <param name="title">グラフタイトル</param>
        /// <param name="yAxisTitle">Y軸タイトル</param>
        /// <param name="fileBaseName">ファイル基本名</param>
        public GraphCreator(GraphConfig config, string title, string yAxisTitle, string fileBaseName)
        {
            _config = config;
            _title = title;
            _yAxisTitle = yAxisTitle;
            _fileBaseName = fileBaseName;
        }


        /// <summary>
        /// 指定した染色体のグラフを作成する。
        /// </summary>
        /// <param name="outDir">出力ディレクトリ</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        /// <returns>グラフファイル</returns>
        public GraphFile Create(OutputDirectory outDir, string chrName, GraphData data)
        {
            var plotModel = CreateNoneSeriesPlotModel(chrName);
            AddSeries(plotModel, chrName, data);
            AddAnnotations(plotModel);

            var filePath = outDir.CreateFilePath($"{_fileBaseName}_{chrName}.png");
            return SaveImageFile(filePath, plotModel);
        }

        /// <summary>
        /// Y軸幅を取得する。
        /// </summary>
        /// <returns></returns>
        public YAxisWidth GetYAxisWidth()
        {
            var dummyPath = $"{_fileBaseName}_dummy.png";
            var plotModel = CreateNoneSeriesPlotModel("dummy");
            var graphFile = SaveImageFile(dummyPath, plotModel);
            graphFile.Delete();

            var yAxis = plotModel.DefaultYAxis;
            return new YAxisWidth(_config.YAxisConfig, yAxis.ScreenMin.X);
        }

        /// <summary>
        /// PlotModelにシリーズを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        /// <param name="chrName">染色体名</param>
        /// <param name="data">グラフデータ</param>
        protected abstract void AddSeries(PlotModel plotModel, string chrName, GraphData data);

        /// <summary>
        /// PlotModelにアノテーションを追加する。
        /// </summary>
        /// <param name="plotModel">PlotModel</param>
        protected abstract void AddAnnotations(PlotModel plotModel);

        /// <summary>
        /// シリーズが空のPlotModelを作成する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>PlotModel</returns>
        private PlotModel CreateNoneSeriesPlotModel(string chrName)
        {
            var plotModel = new PlotModel()
            {
                Title = $"{_title} [{chrName}]",
                TitleFontSize = 25,
                Background = ColorPalette.GraphBackgroundColor
            };

            var xAxis = _config.XAxisConfig.CreateAxis();
            plotModel.Axes.Add(xAxis);

            var yAxis = _config.YAxisConfig.CreateAxis(_yAxisTitle);
            plotModel.Axes.Add(yAxis);

            return plotModel;
        }

        /// <summary>
        /// グラフ画像を保存する。
        /// </summary>
        /// <param name="filePath">グラフ画像ファイルのPath</param>
        /// <param name="plotModel">PlotModel</param>
        /// <returns>グラフファイル</returns>
        private GraphFile SaveImageFile(string filePath, PlotModel plotModel)
        {
            PngExporter.Export(plotModel, filePath, _config.FigureWidth.Value, _config.FigureHeight.Value);

            return new GraphFile(filePath);
        }
    }
}
