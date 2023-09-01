using PolyploidQtlSeqCore.QtlAnalysis.Distribution;
using PolyploidQtlSeqCore.QtlAnalysis.IO;
using PolyploidQtlSeqCore.QtlAnalysis.OxyGraph;
using PolyploidQtlSeqCore.QtlAnalysis.QtlSeqTargetFilter;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオ設定
    /// </summary>
    internal class QtlAnalysisScenarioSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary;

        /// <summary>
        /// コンストラクタ
        /// </summary>
        static QtlAnalysisScenarioSettings()
        {
            var toLongNameDictionary = new Dictionary<string, string>()
            {
                [OutputDirectory.SHORT_NAME] = OutputDirectory.LONG_NAME,
                [OutputDirectory.LONG_NAME] = OutputDirectory.LONG_NAME,

                [DisplayAnnotationImpacts.SHORT_NAME] = DisplayAnnotationImpacts.LONG_NAME,
                [DisplayAnnotationImpacts.LONG_NAME] = DisplayAnnotationImpacts.LONG_NAME,

                [ThreadNumber.SHORT_NAME] = ThreadNumber.LONG_NAME,
                [ThreadNumber.LONG_NAME] = ThreadNumber.LONG_NAME
            };

            QtlSeqTargetPolicySettings.AddLongNameKeyValuePair(toLongNameDictionary);
            NoQtlDistributionSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            SlidingWindowAnalysisSettings.AddLongNameKeyValuePair(toLongNameDictionary);
            GraphSettings.AddLongNameKeyValuePair(toLongNameDictionary);

            _toLongNameDictionary = toLongNameDictionary;
        }

        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// QTL解析シナリオ設定を作成する。
        /// </summary>
        /// <param name="optionValues">オプションの値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション辞書</param>
        public QtlAnalysisScenarioSettings(IQtlAnalysisScenarioSettingValue optionValues, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            OutputDir = new OutputDirectory(optionValues.OutputDir, parameterDictionary, userOptionDictionary);
            DisplayAnnotationImpacts = new DisplayAnnotationImpacts(optionValues.DisplayAnnotationImpacts, parameterDictionary, userOptionDictionary);
            ThreadNumber = new ThreadNumber(optionValues.ThreadNumber, parameterDictionary, userOptionDictionary);

            QtlSeqTargetPolicyOption = new QtlSeqTargetPolicySettings(optionValues, parameterDictionary, userOptionDictionary);
            NoQtlDistributionOption = new NoQtlDistributionSettings(optionValues, parameterDictionary, userOptionDictionary);
            SlidingWindowAnalysisOption = new SlidingWindowAnalysisSettings(optionValues, parameterDictionary, userOptionDictionary);
            GraphOption = new GraphSettings(optionValues, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDir { get; }

        /// <summary>
        /// 表示するAnnotation Impact
        /// </summary>
        public DisplayAnnotationImpacts DisplayAnnotationImpacts { get; }

        /// <summary>
        /// スレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// QTL-seq対象ポリシーオプション
        /// </summary>
        public QtlSeqTargetPolicySettings QtlSeqTargetPolicyOption { get; }

        /// <summary>
        /// QTLなし分布オプション
        /// </summary>
        public NoQtlDistributionSettings NoQtlDistributionOption { get; }

        /// <summary>
        /// SlidingWindow解析オプション
        /// </summary>
        public SlidingWindowAnalysisSettings SlidingWindowAnalysisOption { get; }

        /// <summary>
        /// グラフオプション
        /// </summary>
        public GraphSettings GraphOption { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[] { OutputDir.ToParameterFileLine() }
                .Concat(QtlSeqTargetPolicyOption.ToParameterFileLines())
                .Concat(NoQtlDistributionOption.ToParameterFileLines())
                .Concat(SlidingWindowAnalysisOption.ToParameterFileLines())
                .Concat(GraphOption.ToParameterFileLines())
                .Concat(new[]
                {
                    DisplayAnnotationImpacts.ToParameterFileLine(),
                    ThreadNumber.ToParameterFileLine()
                })
                .ToArray();
        }
    }
}
