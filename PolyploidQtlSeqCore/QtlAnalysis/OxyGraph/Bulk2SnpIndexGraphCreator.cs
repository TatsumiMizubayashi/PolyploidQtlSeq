using OxyPlot;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Bulk2 SnpIndexグラフクリエーター
    /// </summary>
    internal class Bulk2SnpIndexGraphCreator : GraphCreator
    {
        private const string TITLE = "Bulk2 SNP-index";
        private const string Y_AXIS_TITLE = "SNP-index";
        private const string FILE_BASENAME = "Bulk2 SNP-index";

        /// <summary>
        /// Bulk2 SnpIndexグラフクリエーターを作成する。
        /// </summary>
        /// <param name="config">設定</param>
        public Bulk2SnpIndexGraphCreator(GraphConfig config)
            : base(config, TITLE, Y_AXIS_TITLE, FILE_BASENAME)
        {
        }

        protected override void AddAnnotations(PlotModel plotModel)
        {
            // アノテーションなし
        }

        protected override void AddSeries(PlotModel plotModel, string chrName, GraphData data)
        {
            var snpIndesSeries = data.CreateBulk2SnpIndexScatterSeries(chrName);
            plotModel.Series.Add(snpIndesSeries);

            var zeroSnpIndexSeries = data.CreateBulk2ZeroSnpIndexScatterSeries(chrName);
            plotModel.Series.Add(zeroSnpIndexSeries);

            var averageSnpIndexSeries = data.CreateBulk2AverageSnpIndexLineSeries(chrName);
            foreach (var series in averageSnpIndexSeries)
            {
                plotModel.Series.Add(series);
            }
        }
    }
}
