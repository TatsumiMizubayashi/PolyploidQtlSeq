using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// QtlSeqグラフクリエーター
    /// </summary>
    internal class QtlSeqGraphCreator
    {
        private readonly GraphOption _option;

        /// <summary>
        /// QtlSeqグラフクリエーターを作成する。
        /// </summary>
        /// <param name="option">グラフオプション</param>
        public QtlSeqGraphCreator(GraphOption option)
        {
            _option = option;
        }

        /// <summary>
        /// グラフを作成する。
        /// </summary>
        /// <param name="outDir">出力ディレクトリ</param>
        /// <param name="allVariants">全変異データ</param>
        /// <param name="allWindows">全Windowデータ</param>
        public void Create(OutputDirectory outDir, SnpIndexVariantWithSlidingWindowQtl[] allVariants, Window[] allWindows)
        {
            var graphData = new GraphData(allVariants, allWindows, _option.XAxisMajorStep);
            var graphCreator = GetGrhaphCreators(graphData.GraphAxes);

            foreach (var chrName in graphData.ChrNames)
            {
                var mergeGraphFilePath = outDir.CreateFilePath($"{chrName}.png");

                var files = graphCreator.Select(x => x.Create(outDir, chrName, graphData)).ToArray();
                var graphFiles = new GraphFiles(files);
                graphFiles.VerticalMerge(mergeGraphFilePath, _option);
                graphFiles.Delete();
            }
        }

        /// <summary>
        /// グラフクリエーター配列を取得する。
        /// </summary>
        /// <param name="graphAxes">グラフ軸</param>
        /// <returns>グラフクリエーター配列</returns>
        private GraphCreator[] GetGrhaphCreators(GraphAxes graphAxes)
        {
            var correctAxies = graphAxes.CorrectYAxis(_option);

            var snpIndexConfig = correctAxies.CreateSnpIndexGraphConfig(_option);
            var deltaSnpIndexConfig = correctAxies.CreateDeltaSnpIndexGraphConfig(_option);
            var scoreConfig = correctAxies.CreateWindowScoreGraphConfig(_option);
            var qtlCountConfig = correctAxies.CreateQtlCountGraphConfig(_option);

            return new GraphCreator[]
            {
                new Bulk1SnpIndexGraphCreator(snpIndexConfig),
                new Bulk2SnpIndexGraphCreator(snpIndexConfig),
                new BulkSnpIndexGraphCreator(snpIndexConfig),
                new DeltaSnpIndexGraphCreator(deltaSnpIndexConfig),
                new WindowScoreGraphCreator(scoreConfig),
                new QtlVariantCountGraphCreator(qtlCountConfig)
            };
        }
    }
}
