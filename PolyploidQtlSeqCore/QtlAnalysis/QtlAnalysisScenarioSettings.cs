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
        /// <param name="settingValue">設定値</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定オプション辞書</param>
        public QtlAnalysisScenarioSettings(IQtlAnalysisScenarioSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            OutputDir = new OutputDirectory(settingValue.OutputDir, parameterDictionary, userOptionDictionary);
            DisplayAnnotationImpacts = new DisplayAnnotationImpacts(settingValue.DisplayAnnotationImpacts, parameterDictionary, userOptionDictionary);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber, parameterDictionary, userOptionDictionary);

            QtlSeqTargetPolicySettings = new QtlSeqTargetPolicySettings(settingValue, parameterDictionary, userOptionDictionary);
            NoQtlDistributionSettings = new NoQtlDistributionSettings(settingValue, parameterDictionary, userOptionDictionary);
            SlidingWindowAnalysisSettings = new SlidingWindowAnalysisSettings(settingValue, parameterDictionary, userOptionDictionary);
            GraphSettings = new GraphSettings(settingValue, parameterDictionary, userOptionDictionary);
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
        /// QTL-seq対象ポリシー設定
        /// </summary>
        public QtlSeqTargetPolicySettings QtlSeqTargetPolicySettings { get; }

        /// <summary>
        /// QTLなし分布作成設定
        /// </summary>
        public NoQtlDistributionSettings NoQtlDistributionSettings { get; }

        /// <summary>
        /// SlidingWindow解析設定
        /// </summary>
        public SlidingWindowAnalysisSettings SlidingWindowAnalysisSettings { get; }

        /// <summary>
        /// グラフ設定
        /// </summary>
        public GraphSettings GraphSettings { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[] { OutputDir.ToParameterFileLine() }
                .Concat(QtlSeqTargetPolicySettings.ToParameterFileLines())
                .Concat(NoQtlDistributionSettings.ToParameterFileLines())
                .Concat(SlidingWindowAnalysisSettings.ToParameterFileLines())
                .Concat(GraphSettings.ToParameterFileLines())
                .Concat(new[]
                {
                    DisplayAnnotationImpacts.ToParameterFileLine(),
                    ThreadNumber.ToParameterFileLine()
                })
                .ToArray();
        }
    }
}
