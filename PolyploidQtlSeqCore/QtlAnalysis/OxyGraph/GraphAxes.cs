using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Http.Headers;
using System.Text;
using System.Threading.Tasks;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ軸
    /// </summary>
    internal class GraphAxes
    {
        private readonly XAxisConfig _xConfig;
        private readonly YAxisConfig _bulkYConfig;
        private readonly YAxisConfig _deltaYConfig;
        private readonly YAxisConfig _scoreYConfig;
        private readonly YAxisConfig _countYConfig;

        /// <summary>
        /// グラフ軸インスタンスを作成する。
        /// </summary>
        /// <param name="x">X軸範囲</param>
        /// <param name="snpIndex">SnpIndex Y軸設定</param>
        /// <param name="delta">ΔSnpIndex Y軸設定</param>
        /// <param name="score">WindowScore Y軸設定</param>
        /// <param name="count">QTL数 Y軸設定</param>
        public GraphAxes(XAxisConfig x, YAxisConfig snpIndex, YAxisConfig delta, YAxisConfig score, YAxisConfig count)
        {
            _xConfig = x;
            _bulkYConfig = snpIndex;
            _deltaYConfig = delta;
            _scoreYConfig = score;
            _countYConfig = count;
        }

        /// <summary>
        /// Y軸補正を行う。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>補正済みGraphAxes</returns>
        public GraphAxes CorrectYAxis(GraphOption option)
        {
            var snpIndexGraphConfig = CreateSnpIndexGraphConfig(option);
            var deltaSnpIndexGraphConfig = CreateDeltaSnpIndexGraphConfig(option);
            var scoreGraphConfig = CreateWindowScoreGraphConfig(option);
            var qtlCountGraphConfig = CreateQtlCountGraphConfig(option);

            var creators = new GraphCreator[]
            {
                new BulkSnpIndexGraphCreator(snpIndexGraphConfig),
                new DeltaSnpIndexGraphCreator(deltaSnpIndexGraphConfig),
                new WindowScoreGraphCreator(scoreGraphConfig),
                new QtlVariantCountGraphCreator(qtlCountGraphConfig)
            };

            var yAxisWidths = creators.Select(x => x.GetYAxisWidth()).ToArray();
            var maxWidth = yAxisWidths.MaxBy(x => x.Width) ?? throw new InvalidOperationException("maxWidthがnull");
            var correctYConfigs = yAxisWidths.Select(x => x.Correct(maxWidth)).ToArray();

            return new GraphAxes(
                _xConfig,
                correctYConfigs[0],
                correctYConfigs[1],
                correctYConfigs[2],
                correctYConfigs[3]);
        }

        /// <summary>
        /// Bulk SnpIndexグラフ設定を作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>Bulk SnpIndexグラフ設定</returns>
        public GraphConfig CreateSnpIndexGraphConfig(GraphOption option)
        {
            return new GraphConfig(_xConfig, _bulkYConfig, option);
        }

        /// <summary>
        /// ΔSnpIndexグラフ設定を作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>ΔSnpIndexグラフ設定</returns>
        public GraphConfig CreateDeltaSnpIndexGraphConfig(GraphOption option)
        {
            return new GraphConfig(_xConfig, _deltaYConfig, option);
        }


        /// <summary>
        /// WindowScoreグラフ設定を作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>WindowScoreグラフ設定</returns>
        public GraphConfig CreateWindowScoreGraphConfig(GraphOption option)
        {
            return new GraphConfig(_xConfig, _scoreYConfig, option);
        }

        /// <summary>
        /// QTL数グラフ設定を作成する。
        /// </summary>
        /// <param name="option">オプション</param>
        /// <returns>QTL数グラフ設定</returns>
        public GraphConfig CreateQtlCountGraphConfig(GraphOption option)
        {
            return new GraphConfig(_xConfig, _countYConfig, option);
        }
    }
}
