using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.IO
{
    /// <summary>
    /// 表示するアノテーションImpact
    /// </summary>
    public class DisplayAnnotationImpacts
    {
        private static readonly char[] _splitter = new[] { ',' };

        /// <summary>
        /// DisplayImpactsの規定値
        /// </summary>
        public const string DEFAULT = "HIGH,MODERATE";

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "di";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "displayImpacts";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Annotation Impact to be included in the SNP-index file. Separate multiple items with commas.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -di option must be one of HIGH, MODERATE, LOW, or MODIFIER, separated by commas.";

        /// <summary>
        /// 表示するアノテーションImpactを作成する。
        /// </summary>
        /// <param name="impacts">表示するImpact</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public DisplayAnnotationImpacts(string impacts, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Impacts = OptionValue.GetValue(LONG_NAME, impacts, parameterDictionary, userOptionDictionary);

            if (string.IsNullOrEmpty(Impacts)) throw new ArgumentException(null, nameof(impacts));

            var items = Impacts.Split(_splitter);
            var flag = Impact.None;
            foreach (var item in items)
            {
                if (!Enum.TryParse<Impact>(item, out var impact))
                    throw new ArgumentException($"{item} is an invalid value.");

                flag |= impact;
            }

            DisplayFlag = flag;
        }

        /// <summary>
        /// 表示するImapctsのテキストを取得する。
        /// </summary>
        private string Impacts { get; }

        /// <summary>
        /// 表示するImpactフラグを取得する。
        /// Noneの場合はアノテーション情報カラムを表示させない。
        /// </summary>
        internal Impact DisplayFlag { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Impacts}";
        }
    }
}
