using PolyploidQtlSeqCore.QtlAnalysis;
using PolyploidQtlSeqCore.Share;

namespace PolyploidQtlSeqCore.VariantCall
{
    /// <summary>
    /// Bcftools 変異検出オプション
    /// </summary>
    internal class BcfToolsVariantCallSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [MinmumBaseQuality.SHORT_NAME] = MinmumBaseQuality.LONG_NAME,
            [MinmumBaseQuality.LONG_NAME] = MinmumBaseQuality.LONG_NAME,

            [MinmumMappingQuality.SHORT_NAME] = MinmumMappingQuality.LONG_NAME,
            [MinmumMappingQuality.LONG_NAME] = MinmumMappingQuality.LONG_NAME,

            [AdjustMappingQuality.SHORT_NAME] = AdjustMappingQuality.LONG_NAME,
            [AdjustMappingQuality.LONG_NAME] = AdjustMappingQuality.LONG_NAME
        };


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
        /// Bcftools変異検出オプションを作成する。
        /// </summary>
        /// <param name="settingValue">bcftools変異検出設定値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public BcfToolsVariantCallSettings(IBcftoolsVariantCallSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            ReferenceSequence = new ReferenceSequence(settingValue.ReferenceSequence);
            MinmumBaseQuality = new MinmumBaseQuality(settingValue.MinBq, parameterDictionary, userOptionDictionary);
            MinmumMappingQuality = new MinmumMappingQuality(settingValue.MinMq, parameterDictionary, userOptionDictionary);
            AdjustMappingQuality = new AdjustMappingQuality(settingValue.AdjustMq, parameterDictionary, userOptionDictionary);
        }

        /// <summary>
        /// Bcftools変異検出オプションを作成する。
        /// </summary>
        /// <param name="settingValue">bcftools変異検出設定値</param>
        public BcfToolsVariantCallSettings(IBcftoolsVariantCallSettingValue settingValue)
        {
            ReferenceSequence = new ReferenceSequence(settingValue.ReferenceSequence);
            MinmumBaseQuality = new MinmumBaseQuality(settingValue.MinBq);
            MinmumMappingQuality = new MinmumMappingQuality(settingValue.MinMq);
            AdjustMappingQuality = new AdjustMappingQuality(settingValue.AdjustMq);
            OutputDirectory = new OutputDirectory(settingValue.OutputDir);
            ThreadNumber = new ThreadNumber(settingValue.ThreadNumber);
        }

        /// <summary>
        /// リファレンスシークエンス
        /// </summary>
        public ReferenceSequence ReferenceSequence { get; }

        /// <summary>
        /// Base Quality最低値
        /// </summary>
        public MinmumBaseQuality MinmumBaseQuality { get; }

        /// <summary>
        /// Mapping Quality最低値
        /// </summary>
        public MinmumMappingQuality MinmumMappingQuality { get; }

        /// <summary>
        /// Adjust Mapping Quality
        /// </summary>
        public AdjustMappingQuality AdjustMappingQuality { get; }

        /// <summary>
        /// 出力ディレクトリ
        /// </summary>
        public OutputDirectory OutputDirectory { get; }

        /// <summary>
        /// スレッド数
        /// </summary>
        public ThreadNumber ThreadNumber { get; }

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                MinmumBaseQuality.ToParameterFileLine(),
                MinmumMappingQuality.ToParameterFileLine(),
                AdjustMappingQuality.ToParameterFileLine()
            };
        }
    }
}
