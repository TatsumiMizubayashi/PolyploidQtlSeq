using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QualityControl
{
    /// <summary>
    /// リードの最低長
    /// </summary>
    public class ReadLengthRequired
    {
        /// <summary>
        /// リード最低長の最小値
        /// </summary>
        [Obsolete("削除予定")]
        public const int MINIMUM = 10;

        /// <summary>
        /// リード最低長の最大値
        /// </summary>
        [Obsolete("削除予定")]
        public const int MAXIMUM = 150;

        /// <summary>
        /// リード最低長の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 50;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "l";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "lengthRequired";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Minimum lead length after trimming. length_required in fastp.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -l option must be an integer greater than 10 and less than or equal to 150.";


        /// <summary>
        /// リードの最低長を作成する。
        /// </summary>
        /// <param name="length">リード最低長</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        [Obsolete("削除予定")]
        public ReadLengthRequired(int length, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, length, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// リードの最低長を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// Fastpの引数に変換する。
        /// </summary>
        /// <returns>Fastp引数</returns>
        internal string ToFastpArg()
        {
            return $"-l {Value}";
        }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        [Obsolete("削除予定")]
        public string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }

    }
}
